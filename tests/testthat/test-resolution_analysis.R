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

test_that("the same population is sampled at every resolution: one table per resolution", {
  dir <- file.path(tempdir(), "res_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 6, n_snps = 200, n_points = 6), dir)
  reg <- default_ppi_registry()[c("RN", "CV_t")]
  tabs <- suppressMessages(resolution_tables(dir, c(3, 6), reg))
  expect_equal(names(tabs), c("n3", "n6"))
  expect_equal(tabs$n3$Genotype, tabs$n6$Genotype)
  expect_equal(unique(tabs$n3$N_points), 3)
  expect_true(all(c("Shape", "Genotype", "N_points", "RN", "CV_t", summary_stat_names) %in% names(tabs$n6)))
  # at the population's own resolution the PPIs equal those of the population's data
  pop <- read_population(dir)
  rn <- suppressMessages(compute_ppis(pop$phenotypes, pop$environments, reg))$RN
  expect_equal(tabs$n6$RN, rn$Score[match(paste(tabs$n6$Shape, tabs$n6$Genotype), paste(rn$Type, rn$Genotype_ID))])
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

test_that("score table of a population at its own resolution: stored PPIs and summary statistics", {
  dir <- file.path(tempdir(), "score_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200), dir)
  expect_error(population_score_table(dir), "compute_population_ppis")
  suppressMessages(compute_population_ppis(dir, registry = default_ppi_registry()[c("RN", "CV_t")]))
  tab <- population_score_table(dir)
  expect_equal(names(tab), c("Shape", "Genotype", "RN", "CV_t", summary_stat_names))
  expect_equal(nrow(tab), 20L)
  expect_equal(unique(tab$Shape), c("linear", "gaussian", "sinusoidal", "wave"))
})
