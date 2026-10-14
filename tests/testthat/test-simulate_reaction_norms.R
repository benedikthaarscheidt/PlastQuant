small <- function(h2 = 0.4, ...) simulate_reaction_norms(sim_config(h2 = h2, Mc = 5, n_per_shape = 300, n_snps = 2000, seed = 11, ...))

test_that("the result is a valid pipeline input with mother norms and replicates", {
  sim <- small()
  expect_silent(check_pipeline_input(sim))
  d <- sim$phenotypes$gaussian
  expect_equal(names(d), c("Genotype", "Replicate", "Environment", "Trait"))
  expect_equal(nrow(d), 300L * 3L * 6L)
  expect_equal(sim$environments, 1 + 9 * (0:5) / 5)
  expect_equal(unique(d$Genotype)[1], "Genotype_301")
  expect_equal(nrow(sim$mother_curves$wave), 300L * 6L)
  expect_equal(nrow(sim$truth), 19L * 5L)
})

test_that("parameter values stay inside the study space", {
  sim <- small()
  for (s in sim$config$shapes) for (theta in names(sim$config$study_space[[s]])) {
    sp <- sim$config$study_space[[s]][[theta]]
    v <- sim$plant_parameters[[s]][[theta]]
    expect_true(all(v >= sp$min & v <= sp$max), info = paste(s, theta))
  }
})

test_that("common random numbers: only h2 changes between scenarios", {
  a <- small(h2 = 0.2); b <- small(h2 = 0.6)
  expect_identical(a$snp_matrix, b$snp_matrix)
  expect_identical(a$truth, b$truth)
  expect_false(isTRUE(all.equal(a$phenotypes$linear$Trait, b$phenotypes$linear$Trait)))
})

test_that("nothing is negative", {
  sim <- small()
  expect_true(all(vapply(sim$phenotypes, function(d) all(d$Trait > 0), logical(1))))
})

test_that("the reaction norms are generated once and can be sampled at any number of points", {
  cfg <- sim_config(h2 = 0.4, Mc = 2, n_per_shape = 10, n_snps = 300, seed = 5)
  norms <- simulate_norms(cfg)
  expect_null(norms$phenotypes)
  expect_true(all(c("snp_matrix", "truth", "mother_parameters", "plant_parameters", "latent", "config") %in% names(norms)))
  # sampling at the population's own resolution gives exactly the simulated population
  expect_equal(measure_norms(norms, cfg$n_points)$phenotypes, simulate_reaction_norms(cfg)$phenotypes)
  m11 <- measure_norms(norms, 11)
  expect_equal(m11$environments, 1 + 9 * (0:10) / 10)
  expect_equal(nrow(m11$phenotypes$wave), 10L * 3L * 11L)
  # same plants: the noise-free curves agree at the points both designs share (t = 0, 0.2, ..., 1)
  f6 <- measure_norms(norms, 6, error = FALSE)$phenotypes$wave
  f11 <- measure_norms(norms, 11, error = FALSE)$phenotypes$wave
  shared <- f11[f11$Environment %in% f6$Environment, ]
  expect_equal(shared$Trait, f6$Trait[match(paste(shared$Genotype, shared$Replicate, shared$Environment),
                                            paste(f6$Genotype, f6$Replicate, f6$Environment))])
  # the measurement error is the requested share for every design
  expect_true(all(m11$measurement_sd > 0))
  expect_false(isTRUE(all.equal(m11$measurement_sd, measure_norms(norms, 3)$measurement_sd)))
})
