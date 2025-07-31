library(testthat)
library(PlastQuant)  # replace with your actual package name

test_that("returns NA for fewer than 2 values", {
  expect_true(is.na(calculate_D_slope(numeric(0))))
  expect_true(is.na(calculate_D_slope(42)))
})

test_that("errors on missing values in input", {
  expect_error(
    calculate_D_slope(c(1, NA, 3)),
    "missing values"
  )
})

test_that("errors when fractions are out of [0,1]", {
  expect_error(
    calculate_D_slope(1:5, lower_fraction = -0.1),
    "must be between 0 and 1"
  )
  expect_error(
    calculate_D_slope(1:5, upper_fraction = 1.5),
    "must be between 0 and 1"
  )
})

test_that("errors when lower_fraction > upper_fraction", {
  expect_error(
    calculate_D_slope(1:5, lower_fraction = 0.8, upper_fraction = 0.2),
    "`lower_fraction` must be less than `upper_fraction`."
  )
})

test_that("correct D slope with default fractions on sorted data", {
  expect_equal(calculate_D_slope(1:5), 4)
})

test_that("correct D slope with default fractions on unsorted data", {
  vec <- seq(1:10)
  expect_equal(calculate_D_slope(vec), 8)
})

test_that("correct D slope with custom fractions", {
  trait <- 1:5
  expect_equal(calculate_D_slope(trait, lower_fraction = 0.4, upper_fraction = 0.6), 4.5 - 1.5)
})

test_that("handles cases where lower or upper segment has length >1", {
  trait <- 1:10
  expect_equal(calculate_D_slope(trait, lower_fraction = 0.3, upper_fraction = 0.7), 9 - 2)
})
