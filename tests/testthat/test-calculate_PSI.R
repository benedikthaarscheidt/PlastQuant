library(testthat)
library(PlastQuant)

test_that("non-numeric trait_values errors", {
  expect_error(calculate_PSI("a"), "`trait_values` must be a numeric vector")
})

test_that("too few non-NA trait_values errors", {
  expect_error(calculate_PSI(c(NA, 5)), "At least two non-NA trait values")
})

test_that("mismatched env_values length errors", {
  expect_error(calculate_PSI(c(1,2,3), env_values = c(1,2)),
               "must have the same length")
})

test_that("non-numeric env_values errors", {
  expect_error(calculate_PSI(c(1,2,3), env_values = c("a","b","c")),
               "Cannot coerce character `env_values` to numeric")
})

test_that("zero env range returns NA with warning", {
  vec <- c(1,2,3)
  expect_warning(res <- calculate_PSI(vec, env_values = c(5,5,5)),
                 "zero or near-zero range")
  expect_true(is.na(res))
})

test_that("PSI correct for simple positive slope", {
  trait <- c(1, 2, 3, 4)
  env   <- c(1, 2, 3, 4)
  # slope = 1, stability = 1/(1+1) = 0.5
  expect_equal(calculate_PSI(trait, env), 0.5)
})

test_that("PSI correct for negative slope", {
  trait <- c(4,3,2,1)
  env   <- c(1,2,3,4)
  # slope = -1, stability = 1/(1+1) = 0.5
  expect_equal(calculate_PSI(trait, env), 0.5)
})

test_that("PSI handles missing trait_values", {
  trait <- c(NA, 2, NA, 4)
  env   <- c(1, 2, 3, 4)
  # drop NAs => trait=c(2,4), env=c(2,4) slope=1 => 0.5
  expect_equal(calculate_PSI(trait, env), 0.5)
})
