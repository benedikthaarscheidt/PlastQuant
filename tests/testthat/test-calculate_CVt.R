library(testthat)
library(PlastQuant)

test_that("calculate_CVt returns NA for vectors of length < 2", {
  expect_true(is.na(calculate_CVt(numeric(0))))
  expect_true(is.na(calculate_CVt(42)))
})

test_that("calculate_CVt errors on any NA in the input", {
  expect_error(
    calculate_CVt(c(1, NA, 3)),
    "missing values"
  )
})

test_that("calculate_CVt computes sd(x)/mean(x) for normal numeric vectors", {
  x <- c(2, 4, 6)
  # For this x, mean = 4, sd = 2, so CVt = 2/4 = 0.4
  expect_equal(calculate_CVt(x), 0.5)
  # also test with non‑integer values
  y <- c(10.2, 9.8, 10.0)
  expected <- sd(y) / mean(y)
  expect_equal(calculate_CVt(y), expected)
})

test_that("calculate_CVt handles zero mean appropriately", {
  z <- c(-1, 0, 1) * 10  # mean = 0, sd > 0
  expect_true(is.na(calculate_CVt(z)))
})
