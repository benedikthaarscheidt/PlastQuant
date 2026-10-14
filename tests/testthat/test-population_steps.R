fresh_dir <- function(name) {
  d <- file.path(tempdir(), name)
  unlink(d, recursive = TRUE)
  d
}

test_that("simulation, PPIs and GWAS are separate steps on the population folder", {
  dir <- fresh_dir("steps")
  simulate_population(sim_config(h2 = 0.6, Mc = 5, n_per_shape = 40, n_snps = 400, seed = 3), dir)
  expect_true(all(file.exists(file.path(dir, c("data/traits/linear.csv", "data/environments.csv", "data/snp_matrix.rds",
                                               "data/truth.csv", "data/config.rds", "data/mother_parameters.csv",
                                               "variance/variance_latent.csv", "variance/variance_trait.csv",
                                               "variance/measurement.csv")))))
  expect_false(dir.exists(file.path(dir, "ppis")))
  expect_error(run_population_gwas(dir), "compute_population_ppis")

  reg <- default_ppi_registry()[c("CV_t", "RN", "FW")]
  expect_equal(suppressMessages(compute_population_ppis(dir, registry = reg)), dir)
  expect_true(file.exists(file.path(dir, "ppis", "RN.csv")))
  expect_false(dir.exists(file.path(dir, "gwas")))

  expect_equal(suppressMessages(run_population_gwas(dir)), dir)
  power <- read.csv(file.path(dir, "gwas", "gwas_power.csv"))
  expect_setequal(unique(power$Parameter[power$Shape == "linear"]), c("L", "a", "V", "False_Positives"))
  expect_setequal(unique(read.csv(file.path(dir, "gwas", "gwas_results.csv"))$Score), c("CV_t", "RN", "FW"))
})

test_that("recomputing the PPIs overwrites the earlier ones", {
  dir <- fresh_dir("steps_rerun")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 8, n_snps = 200), dir)
  suppressMessages(compute_population_ppis(dir, registry = default_ppi_registry()[c("CV_t", "RN")]))
  suppressMessages(compute_population_ppis(dir, registry = default_ppi_registry()["CV_t"]))
  expect_equal(sort(list.files(file.path(dir, "ppis"))), c("CV_t.csv", "ppi_scores.rds"))
})

test_that("own data without SNPs: PPIs only; the GWAS step names what is missing", {
  ph <- data.frame(Genotype = rep(c("var_A", "var_B", "var_C"), each = 2), Replicate = 1L,
                   Environment = rep(c(5, 9), 3), Trait = c(1, 2, 1, 3, 2, 2))
  dir <- write_population(pipeline_input(list(biomass = ph), c(5, 9)), fresh_dir("own"))
  suppressMessages(compute_population_ppis(dir, registry = default_ppi_registry()["CV_t"]))
  expect_true(file.exists(file.path(dir, "ppis", "CV_t.csv")))
  expect_error(run_population_gwas(dir), "SNP matrix")
})

test_that("a kinship matrix in data/ is used for the GWAS", {
  dir <- fresh_dir("own_kin")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 10, n_snps = 200), dir)
  pop <- read_population(dir)
  K <- kinship_matrix(pop$snp_matrix)
  K[] <- 0; diag(K) <- 1
  write.csv(data.frame(Genotype = rownames(K), K, check.names = FALSE), file.path(dir, "data", "kinship.csv"), row.names = FALSE)
  reg <- default_ppi_registry()["CV_t"]
  suppressMessages(compute_population_ppis(dir, registry = reg))
  suppressMessages(run_population_gwas(dir))
  with_file <- read.csv(file.path(dir, "gwas", "gwas_results.csv"))
  traits <- readRDS(file.path(dir, "ppis", "ppi_scores.rds"))
  expect_equal(with_file$P_Value, run_gwas(traits, pop$snp_matrix, K)$P_Value, tolerance = 1e-8)
})

test_that("the three steps of a scenario, with one GWAS row per replicate", {
  dir <- fresh_dir("reps")
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 30, n_snps = 300, seed = 4)
  suppressMessages({
    simulate_population(cfg, dir)
    compute_population_ppis(dir, registry = default_ppi_registry()["CV_t"], replicate_handling = "replicates")
    run_population_gwas(dir)
  })
  sc <- readRDS(file.path(dir, "ppis", "ppi_scores.rds"))$CV_t
  expect_equal(nrow(sc), 4L * 30L * 3L)
  expect_true(file.exists(file.path(dir, "gwas", "gwas_power.csv")))
})
