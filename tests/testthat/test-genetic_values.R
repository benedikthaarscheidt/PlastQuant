test_that("genetic value has variance 1 and SNP share h2", {
  set.seed(1)
  S <- rnorm(1e5); B <- rnorm(1e5)
  ell <- genetic_value(S, B, 0.3)
  expect_equal(var(ell), 1, tolerance = 0.02)
  expect_equal(var(sqrt(0.3) * S) / var(ell), 0.3, tolerance = 0.02)
})

test_that("repeatability follows Beta(18, 2)", {
  set.seed(2)
  rho <- repeatability_from_latent(rnorm(20000), c(18, 2))
  expect_equal(mean(rho), 0.9, tolerance = 0.005)
  expect_gt(ks.test(rho, "pbeta", 18, 2)$p.value, 0.001)
})

test_that("replicate spread gives correlation rho between replicates", {
  set.seed(3)
  ell <- rnorm(50000); tau <- replicate_spread(0.8)
  z <- matrix(rnorm(100000), 50000, 2)
  p <- plant_latent_values(ell, tau, z)
  expect_equal(cor(p[, 1], p[, 2]), 0.8, tolerance = 0.01)
})

test_that("extreme latent values stay finite", {
  rho <- repeatability_from_latent(c(-40, -9, 9, 40), c(18, 2))
  expect_true(all(rho > 0 & rho < 1))
  expect_true(all(is.finite(replicate_spread(rho))))
})
