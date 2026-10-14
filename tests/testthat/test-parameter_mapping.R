cfg <- sim_config(h2 = 0.4, Mc = 5)

test_that("mother norms are (log-)uniform on the study space", {
  set.seed(1)
  ell <- rnorm(20000)
  expect_gt(ks.test(latent_to_parameter(ell, "a", "linear", cfg), "punif", 0, 0.7)$p.value, 0.001)
  expect_gt(ks.test(log(latent_to_parameter(ell, "w", "gaussian", cfg)), "punif", log(0.05), 0)$p.value, 0.001)
  expect_true(all(range(latent_to_parameter(ell, "k", "wave", cfg)) >= 1))
})

test_that("L is log-normal with sigma_L", {
  set.seed(2)
  expect_equal(sd(log(latent_to_parameter(rnorm(20000), "L", "linear", cfg))), 0.1, tolerance = 0.005)
})

test_that("all maps are increasing and keep matrix shape", {
  ell <- seq(-3, 3, by = 0.1)
  for (x in list(c("L", "gaussian"), c("a", "gaussian"), c("c", "gaussian"), c("w", "gaussian"),
                 c("k", "sinusoidal"), c("psi", "sinusoidal"))) {
    expect_true(all(diff(latent_to_parameter(ell, x[1], x[2], cfg)) > 0), info = x[1])
  }
  expect_equal(dim(latent_to_parameter(matrix(0, 2, 3), "c", "gaussian", cfg)), c(2L, 3L))
})

test_that("extreme values stay inside the range; unknown parameter is an error", {
  psi <- latent_to_parameter(c(-40, 40), "psi", "sinusoidal", cfg)
  expect_true(all(psi > 0 & psi < pi))
  expect_error(latent_to_parameter(0, "w", "linear", cfg), "no study space for parameter w of shape linear")
})
