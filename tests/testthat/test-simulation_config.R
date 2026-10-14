test_that("shapes and their parameters match the model", {
  cfg <- sim_config(h2 = 0.4, Mc = 5)
  expect_equal(cfg$shapes, c("linear", "gaussian", "sinusoidal", "wave"))
  expect_equal(reaction_norm_parameters(cfg, "linear"), c("L", "a"))
  expect_equal(reaction_norm_parameters(cfg, "gaussian"), c("L", "a", "c", "w"))
  expect_equal(reaction_norm_parameters(cfg, "sinusoidal"), c("L", "a", "k", "psi"))
  expect_equal(reaction_norm_parameters(cfg, "wave"), c("L", "a", "k", "c", "w"))
})

test_that("study space matches the table of the model", {
  sp <- default_study_space()
  expect_equal(sp$linear$a, list(type = "uniform", min = 0, max = 0.7))
  expect_equal(sp$gaussian$w, list(type = "loguniform", min = 0.05, max = 1))
  expect_equal(sp$sinusoidal$k, list(type = "loguniform", min = 0.12, max = 2))
  expect_equal(sp$wave$k, list(type = "loguniform", min = 1, max = 2))
  expect_equal(sp$sinusoidal$psi, list(type = "uniform", min = 0, max = pi))
})

test_that("defaults of the variance components", {
  cfg <- sim_config(h2 = 0.2, Mc = 10)
  expect_equal(cfg$sigma_L, 0.1)
  expect_equal(cfg$repeatability_beta, c(18, 2))
  expect_equal(cfg$error_share, 0.05)
})

test_that("invalid settings fail with a clear message", {
  expect_error(sim_config(h2 = 1.2, Mc = 5), "h2 must be in \\[0, 1\\]")
  expect_error(sim_config(h2 = 0.4, Mc = 5, n_points = 1), "n_points must be at least 2")
  expect_error(sim_config(h2 = 0.4, Mc = 5, n_rep = 0), "n_rep must be at least 1")
  expect_error(sim_config(h2 = 0.4, Mc = 5, error_share = 1), "error_share")
})
