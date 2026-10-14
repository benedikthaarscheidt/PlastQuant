test_that("a small scenario runs end to end and writes all results", {
  out <- file.path(tempdir(), "pipeline_test")
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 40, n_snps = 400, seed = 3)
  reg <- default_ppi_registry()[c("CV_t", "RN", "FW")]
  run_dir <- suppressMessages(run_simulation_scenario(cfg, out, registry = reg))
  files <- list.files(run_dir)
  for (f in c("phenotypes.rds", "snp_matrix.rds", "ppi_scores.rds", "truth.csv", "gwas_results.csv",
              "gwas_power.csv", "variance_latent.csv", "variance_trait.csv", "measurement.csv",
              "mother_parameters.csv", "mother_curves.rds", "config.rds")) expect_true(f %in% files, info = f)
  power <- read.csv(file.path(run_dir, "gwas_power.csv"))
  expect_setequal(unique(power$Parameter[power$Shape == "linear"]), c("L", "a", "V", "False_Positives"))
})

test_that("own data without SNPs runs the PPIs only", {
  ph <- data.frame(Genotype = rep(c("var_A", "var_B", "var_C"), each = 2), Replicate = 1L,
                   Environment = rep(c(5, 9), 3), Trait = c(1, 2, 1, 3, 2, 2))
  x <- pipeline_input(list(biomass = ph), c(5, 9))
  out <- suppressMessages(run_pipeline(x, file.path(tempdir(), "own"), registry = default_ppi_registry()["CV_t"], gwas = FALSE))
  expect_true(file.exists(file.path(out, "ppi_scores.rds")))
  expect_error(run_pipeline(x, tempdir()), "gwas = TRUE needs a snp_matrix")
})

test_that("a scenario runs with one GWAS row per replicate", {
  cfg <- sim_config(h2 = 0.6, Mc = 5, n_per_shape = 30, n_snps = 300, seed = 4)
  run_dir <- suppressMessages(run_simulation_scenario(cfg, file.path(tempdir(), "reps"),
                                                      registry = default_ppi_registry()["CV_t"],
                                                      replicate_handling = "replicates"))
  sc <- readRDS(file.path(run_dir, "ppi_scores.rds"))$CV_t
  expect_equal(nrow(sc), 4L * 30L * 3L)
  expect_true(file.exists(file.path(run_dir, "gwas_power.csv")))
})
