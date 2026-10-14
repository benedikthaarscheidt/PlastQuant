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
