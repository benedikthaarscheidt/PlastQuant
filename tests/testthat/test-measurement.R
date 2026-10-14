test_that("measurement error is the requested share of the trait variance", {
  set.seed(1)
  Y <- matrix(rnorm(6000, 1, 0.2), 1000, 6)
  s <- measurement_error_sd(Y, 0.05)
  expect_equal(s^2 / (s^2 + mean(apply(Y, 2, var))), 0.05)
  expect_equal(add_measurement_error(Y, s, matrix(1, 1000, 6)), Y + s)
})
