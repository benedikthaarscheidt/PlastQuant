fresh_dir <- function(name) {
  d <- file.path(tempdir(), name)
  unlink(d, recursive = TRUE)
  d
}

test_that("a simulated population is saved and all results are written into its folder", {
  dir <- fresh_dir("pipeline_test")
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 40, n_snps = 400, seed = 3)
  reg <- default_ppi_registry()[c("CV_t", "RN", "FW")]
  expect_equal(suppressMessages(run_simulation_scenario(cfg, dir, registry = reg)), dir)
  for (f in c("data/traits/linear.csv", "data/environments.csv", "data/snp_matrix.rds", "data/truth.csv",
              "data/config.rds", "data/mother_parameters.csv", "ppis/ppi_scores.rds", "ppis/RN.csv",
              "gwas/gwas_results.csv", "gwas/gwas_power.csv", "variance/variance_latent.csv",
              "variance/variance_trait.csv", "variance/measurement.csv")) {
    expect_true(file.exists(file.path(dir, f)), info = f)
  }
  power <- read.csv(file.path(dir, "gwas", "gwas_power.csv"))
  expect_setequal(unique(power$Parameter[power$Shape == "linear"]), c("L", "a", "V", "False_Positives"))
})

test_that("rerunning the pipeline overwrites the results of the population", {
  dir <- fresh_dir("pipeline_rerun")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 8, n_snps = 200), dir)
  suppressMessages(run_pipeline(dir, registry = default_ppi_registry()[c("CV_t", "RN")], gwas = FALSE))
  suppressMessages(run_pipeline(dir, registry = default_ppi_registry()["CV_t"], gwas = FALSE))
  expect_equal(sort(list.files(file.path(dir, "ppis"))), c("CV_t.csv", "ppi_scores.rds"))
})

test_that("own data without SNPs runs the PPIs only", {
  ph <- data.frame(Genotype = rep(c("var_A", "var_B", "var_C"), each = 2), Replicate = 1L,
                   Environment = rep(c(5, 9), 3), Trait = c(1, 2, 1, 3, 2, 2))
  dir <- write_population(pipeline_input(list(biomass = ph), c(5, 9)), fresh_dir("own"))
  suppressMessages(run_pipeline(dir, registry = default_ppi_registry()["CV_t"]))
  expect_true(file.exists(file.path(dir, "ppis", "CV_t.csv")))
  expect_false(dir.exists(file.path(dir, "gwas")))
  expect_error(run_pipeline(dir, gwas = TRUE), "gwas = TRUE needs a SNP matrix")
})

test_that("a kinship matrix in data/ is used for the GWAS", {
  dir <- fresh_dir("own_kin")
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 10, n_snps = 200), dir)
  pop <- read_population(dir)
  K <- kinship_matrix(pop$snp_matrix)
  K[] <- 0; diag(K) <- 1
  write.csv(data.frame(Genotype = rownames(K), K, check.names = FALSE), file.path(dir, "data", "kinship.csv"), row.names = FALSE)
  reg <- default_ppi_registry()["CV_t"]
  suppressMessages(run_pipeline(dir, registry = reg))
  with_file <- read.csv(file.path(dir, "gwas", "gwas_results.csv"))
  traits <- suppressMessages(compute_ppis(pop$phenotypes, pop$environments, reg))
  expected <- run_gwas(traits, pop$snp_matrix, K)
  expect_equal(with_file$P_Value, expected$P_Value, tolerance = 1e-8)
})

test_that("a scenario runs with one GWAS row per replicate", {
  dir <- fresh_dir("reps")
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 30, n_snps = 300, seed = 4)
  suppressMessages(run_simulation_scenario(cfg, dir, registry = default_ppi_registry()["CV_t"],
                                           replicate_handling = "replicates"))
  sc <- readRDS(file.path(dir, "ppis", "ppi_scores.rds"))$CV_t
  expect_equal(nrow(sc), 4L * 30L * 3L)
  expect_true(file.exists(file.path(dir, "gwas", "gwas_power.csv")))
})
