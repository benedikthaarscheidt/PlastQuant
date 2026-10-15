fresh_dir <- function(name) {
  d <- file.path(tempdir(), name)
  unlink(d, recursive = TRUE)
  d
}

test_that("simulation, sampling at a resolution, PPIs and GWAS are separate steps", {
  dir <- fresh_dir("steps")
  simulate_population(sim_config(h2 = 0.6, Mc = 5, n_per_shape = 40, n_snps = 400, seed = 3), dir)
  expect_true(all(file.exists(file.path(dir, c("data/norms.rds", "data/snp_matrix.rds", "data/truth.csv",
                                               "data/config.rds", "data/mother_parameters.csv",
                                               "variance/variance_latent.csv")))))
  expect_length(resolution_dirs(dir), 0)

  res <- sample_resolution(dir, 6)
  expect_equal(res, file.path(dir, "n_6"))
  expect_true(all(file.exists(file.path(res, c("data/traits/linear.csv", "data/environments.csv",
                                               "variance/variance_trait.csv", "variance/measurement.csv")))))
  expect_error(run_population_gwas(res), "compute_population_ppis")

  reg <- default_ppi_registry()[c("CV_t", "RN", "FW")]
  expect_equal(suppressMessages(compute_population_ppis(res, registry = reg)), res)
  expect_true(file.exists(file.path(res, "mean_before", "ppis", "RN.csv")))
  expect_false(dir.exists(file.path(res, "mean_before", "gwas")))

  expect_equal(suppressMessages(run_population_gwas(res)), res)
  power <- read.csv(file.path(res, "mean_before", "gwas", "gwas_power.csv"))
  expect_setequal(unique(power$Parameter[power$Shape == "linear"]), c("L", "a", "V", "False_Positives"))
  expect_setequal(unique(read.csv(file.path(res, "mean_before", "gwas", "gwas_results.csv"))$Score), c("CV_t", "RN", "FW"))
})

test_that("the same population at two resolutions: same plants, separate folders", {
  dir <- fresh_dir("steps_two")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200), dir)
  r6 <- sample_resolution(dir, 6)
  r11 <- sample_resolution(dir, 11)
  expect_equal(basename(resolution_dirs(dir)), c("n_6", "n_11"))
  p6 <- read_population(r6); p11 <- read_population(r11)
  expect_equal(length(p11$environments), 11L)
  expect_equal(unique(p6$phenotypes$wave$Genotype), unique(p11$phenotypes$wave$Genotype))
  expect_identical(p6$snp_matrix, p11$snp_matrix)
  # sampling again gives the same measurement and clears the results of the old one
  suppressMessages(compute_population_ppis(r6, registry = default_ppi_registry()["CV_t"]))
  sample_resolution(dir, 6)
  expect_false(dir.exists(file.path(r6, "mean_before")))
  expect_equal(read_population(r6)$phenotypes, p6$phenotypes)
})

test_that("recomputing the PPIs overwrites the earlier ones", {
  dir <- fresh_dir("steps_rerun")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 8, n_snps = 200), dir)
  res <- sample_resolution(dir, 6)
  suppressMessages(compute_population_ppis(res, registry = default_ppi_registry()[c("CV_t", "RN")]))
  suppressMessages(compute_population_ppis(res, registry = default_ppi_registry()["CV_t"]))
  expect_equal(sort(list.files(file.path(res, "mean_before", "ppis"))), c("CV_t.csv", "ppi_scores.rds"))
})

test_that("own data without SNPs: PPIs only; the GWAS step names what is missing", {
  ph <- data.frame(Genotype = rep(c("var_A", "var_B", "var_C"), each = 2), Replicate = 1L,
                   Environment = rep(c(5, 9), 3), Trait = c(1, 2, 1, 3, 2, 2))
  dir <- write_population(pipeline_input(list(biomass = ph), c(5, 9)), fresh_dir("own"))
  res <- resolution_dirs(dir)
  suppressMessages(compute_population_ppis(res, registry = default_ppi_registry()["CV_t"]))
  expect_true(file.exists(file.path(res, "mean_before", "ppis", "CV_t.csv")))
  expect_error(run_population_gwas(res), "SNP matrix")
  expect_error(sample_resolution(dir, 3), "simulated population")
})

test_that("a kinship matrix in the population's data/ is used for the GWAS", {
  dir <- fresh_dir("own_kin")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 10, n_snps = 200), dir)
  res <- sample_resolution(dir, 6)
  pop <- read_population(res)
  K <- kinship_matrix(pop$snp_matrix)
  K[] <- 0; diag(K) <- 1
  write.csv(data.frame(Genotype = rownames(K), K, check.names = FALSE), file.path(dir, "data", "kinship.csv"), row.names = FALSE)
  reg <- default_ppi_registry()["CV_t"]
  suppressMessages(compute_population_ppis(res, registry = reg))
  suppressMessages(run_population_gwas(res))
  with_file <- read.csv(file.path(res, "mean_before", "gwas", "gwas_results.csv"))
  traits <- readRDS(file.path(res, "mean_before", "ppis", "ppi_scores.rds"))
  expect_equal(with_file$P_Value, run_gwas(traits, pop$snp_matrix, K)$P_Value, tolerance = 1e-8)
})

test_that("the steps of a scenario, with one GWAS row per replicate", {
  dir <- fresh_dir("reps")
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 30, n_snps = 300, seed = 4)
  suppressMessages({
    simulate_population(cfg, dir)
    res <- sample_resolution(dir, 6)
    compute_population_ppis(res, registry = default_ppi_registry()["CV_t"], replicate_handling = "replicates")
    run_population_gwas(res, replicate_handling = "replicates")
  })
  sc <- readRDS(file.path(res, "replicates", "ppis", "ppi_scores.rds"))$CV_t
  expect_equal(nrow(sc), 4L * 30L * 3L)
  expect_true(file.exists(file.path(res, "replicates", "gwas", "gwas_power.csv")))
})

test_that("PPIs from averaged and from single replicates live side by side, each with its own GWAS", {
  dir <- fresh_dir("both_handlings")
  simulate_population(sim_config(h2 = 0.6, Mc = 2, n_per_shape = 10, n_snps = 200), dir)
  res <- sample_resolution(dir, 6)
  reg <- default_ppi_registry()["CV_t"]
  suppressMessages({
    compute_population_ppis(res, registry = reg)
    compute_population_ppis(res, registry = reg, replicate_handling = "replicates")
    run_population_gwas(res)
    run_population_gwas(res, replicate_handling = "replicates")
  })
  avg <- readRDS(file.path(res, "mean_before", "ppis", "ppi_scores.rds"))$CV_t
  per <- readRDS(file.path(res, "replicates", "ppis", "ppi_scores.rds"))$CV_t
  expect_equal(nrow(avg), 40L)
  expect_equal(nrow(per), 120L)
  expect_true(all(file.exists(file.path(res, c("mean_before", "replicates"), "gwas", "gwas_results.csv"))))
  expect_error(run_population_gwas(file.path(dir, "n_6"), replicate_handling = "averaged"), "should be one of")
})
