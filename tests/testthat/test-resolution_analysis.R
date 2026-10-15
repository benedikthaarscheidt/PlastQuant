test_that("summary statistics of a reaction norm, slope on the environment values", {
  M <- rbind(g1 = c(1, 2, 3, 4), g2 = c(4, 4, 4, 4))
  st <- summary_statistics(M, environments = c(0, 2, 4, 6))
  expect_equal(names(st), c("Genotype", summary_stat_names))
  expect_equal(st$slope, c(0.5, 0))
  expect_equal(st$mean, c(2.5, 4))
  expect_equal(st$mean_lower, c(1.5, 4))
  expect_equal(st$mean_upper, c(3.5, 4))
  expect_equal(st$range, c(3, 0))
  expect_equal(st$variance, c(var(1:4), 0))
  expect_equal(st$median, c(2.5, 4))
})

test_that("a sampling design is the old index selection on a grid of 50 points", {
  expect_equal(grid_design(1, 50, 10), c(1, 11, 21, 31, 41, 50))
  expect_equal(grid_design(1, 50, 49), c(1, 50))
  expect_equal(grid_design(25, 50, 20), c(25, 45, 50))
  expect_equal(design_points(c(1, 50), 50), c(0, 1))
  expect_equal(design_points(c(1, 25), 50), c(0, 24 / 49))
})

test_that("all resolutions of one population: sampled, PPIs and summary statistics in their n_* folders", {
  dir <- file.path(tempdir(), "res_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 6, n_snps = 200), dir)
  reg <- default_ppi_registry()[c("RN", "CV_t")]
  res <- suppressMessages(compute_resolutions(dir, c(3, 6), reg))
  expect_equal(basename(res), c("n_3", "n_6"))
  for (r in res) expect_true(all(file.exists(file.path(r, c("data/environments.csv", "mean_before/ppis/ppi_scores.rds",
                                                             "analyses/summary_statistics.csv")))))
  # a resolution that already has its PPIs is not computed again
  stamp <- file.mtime(file.path(res[2], "mean_before", "ppis", "ppi_scores.rds"))
  Sys.sleep(1.1)
  suppressMessages(compute_resolutions(dir, 6, reg))
  expect_equal(file.mtime(file.path(res[2], "mean_before", "ppis", "ppi_scores.rds")), stamp)

  tabs <- resolution_score_tables(dir)
  expect_equal(names(tabs), c("n3", "n6"))
  expect_equal(unique(tabs$n6$N_points), 6)
  expect_equal(tabs$n3$Genotype, tabs$n6$Genotype)
  expect_true(all(c("Shape", "Genotype", "N_points", "RN", "CV_t", summary_stat_names) %in% names(tabs$n6)))
  expect_equal(names(resolution_score_tables(dir, n_points = 6)), "n6")
  # the PPIs are those of the resolution folder
  rn <- readRDS(file.path(res[2], "mean_before", "ppis", "ppi_scores.rds"))$RN
  expect_equal(tabs$n6$RN, rn$Score[match(paste(tabs$n6$Shape, tabs$n6$Genotype), paste(rn$Type, rn$Genotype_ID))])
  # layout of the earlier regression-data tables: genotype, scores, statistics
  legacy <- as_regression_data(tabs$n6)
  expect_equal(names(legacy), c("Genotype", "RN", "CV_t", summary_stat_names))
  expect_equal(legacy$Genotype[1], "linear_Genotype_1")
})

toy_tables <- function(n_set = c(2, 4, 6, 8), n_gen = 40, seed = 1) {
  set.seed(seed)
  base <- rnorm(n_gen)
  setNames(lapply(n_set, function(n) {
    st <- as.data.frame(matrix(rnorm(n_gen * 9), n_gen, dimnames = list(NULL, summary_stat_names)))
    st$mean <- base + rnorm(n_gen, sd = 0.01)
    data.frame(Shape = "linear", Genotype = paste0("G", seq_len(n_gen)), N_points = n,
               depends = 0.5 * n + base + rnorm(n_gen, sd = 0.1),
               constant = base + rnorm(n_gen, sd = 0.1), st)
  }), paste0("n", n_set))
}

test_that("Pearson correlation of every PPI with every summary statistic, per resolution", {
  cors <- correlations_by_resolution(toy_tables())
  expect_equal(names(cors), c("N_points", "PPI", "Stat", "r"))
  expect_equal(nrow(cors), 4L * 2L * 9L)
  expect_true(all(cors$r[cors$PPI == "constant" & cors$Stat == "mean"] > 0.9))
})

test_that("regression of the correlation on the number of points", {
  cors <- data.frame(N_points = rep(c(2, 4, 6, 8), 2), PPI = rep(c("A", "B"), each = 4), Stat = "mean",
                     r = c(0.2, 0.4, 0.6, 0.8 + 1e-3, 0.5, 0.5, NA, NA))
  reg <- correlation_regression(cors)
  expect_equal(reg$slope[reg$PPI == "A"], 0.1, tolerance = 1e-3)
  expect_true(reg$p_value[reg$PPI == "A"] < 0.01)
  expect_true(is.na(reg$slope[reg$PPI == "B"]))     # fewer than 3 resolutions with a value
})

test_that("mixed model: effect of the number of points on a PPI after the summary statistics", {
  mm <- sample_size_mixed_model(toy_tables())
  expect_equal(sort(mm$PPI), c("constant", "depends"))
  expect_equal(mm$estimate[mm$PPI == "depends"], 0.5, tolerance = 0.02)
  expect_equal(mm$class[mm$PPI == "depends"], "sample-size-sensitive")
  expect_equal(mm$class[mm$PPI == "constant"], "robust")
  expect_true(all(c("std.error", "t.value", "p.value", "p.adj", "ci_lower", "ci_upper") %in% names(mm)))
})

test_that("score table of one resolution: stored PPIs and summary statistics", {
  dir <- file.path(tempdir(), "score_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200), dir)
  res <- sample_resolution(dir, 6)
  expect_error(population_score_table(res), "compute_population_ppis")
  suppressMessages(compute_population_ppis(res, registry = default_ppi_registry()[c("RN", "CV_t")]))
  tab <- population_score_table(res)
  expect_equal(names(tab), c("Shape", "Genotype", "RN", "CV_t", summary_stat_names))
  expect_equal(nrow(tab), 20L)
  expect_equal(unique(tab$Shape), c("linear", "gaussian", "sinusoidal", "wave"))
})

test_that("inputs of the score analyses from a resolution folder: PPIs, mean reaction norms, shapes", {
  dir <- file.path(tempdir(), "inputs_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 4, n_snps = 200), dir)
  res <- sample_resolution(dir, 6)
  expect_error(score_analysis_inputs(res), "compute_population_ppis")
  suppressMessages(compute_population_ppis(res, registry = default_ppi_registry()[c("RN", "FW")]))
  inp <- score_analysis_inputs(res)
  expect_equal(inp$genotype_forms, c("linear", "gaussian", "sinusoidal", "wave"))
  expect_equal(names(inp$scores_list), c("RN", "FW"))
  expect_equal(dim(inp$rn_matrix_list$wave), c(4L, 6L))
  # rows of the matrices in the order of the scores of that shape
  sc <- inp$scores_list$RN
  expect_equal(rownames(inp$rn_matrix_list$gaussian), sc$Genotype_ID[sc$Type == "gaussian"])
})
