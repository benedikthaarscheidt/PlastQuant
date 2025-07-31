library(testthat)
library(PlastQuant)

# Length mismatch between trait and env triggers error
test_that("length mismatch triggers error", {
  trait <- c(1,2,3)
  env   <- c(0.1, 0.2)
  expect_error(
    calculate_RTR(trait, env),
    "must have the same length"
  )
})

# Non-numeric inputs trigger error
test_that("non-numeric inputs trigger error", {
  expect_error(
    calculate_RTR("a", 1:3),
    "must be numeric vectors"
  )
  expect_error(
    calculate_RTR(1:3, "a"),
    "must be numeric vectors"
  )
})

# All NA after removing NAs returns NA with warning
test_that("all NA returns NA with warning", {
  trait <- c(NA, NA)
  env   <- c(NA, NA)
  expect_warning(
    res <- calculate_RTR(trait, env),
    "All observations have NA"
  )
  expect_true(is.na(res))
})

# Denominator zero (all trait zero) returns NA with warning
test_that("denominator zero returns NA with warning", {
  trait <- c(0, 0, 0, 0)
  env   <- c(1, 2, 3, 4)
  expect_warning(
    res <- calculate_RTR(trait, env),
    "Maximum absolute trait value is zero"
  )
  expect_true(is.na(res))
})


# Quantile-based thresholds work correctly
test_that("quantile thresholds compute correct RTR", {
  trait <- 1:100
  env   <- 1:100
  # lowest 10% = env <= 10, highest 10% = env >= 91
  expected <- (mean(trait[91:100]) - mean(trait[1:10])) / max(trait)
  res <- calculate_RTR(trait, env, env_low = 0.1, env_high = 0.1)
  expect_equal(res, expected)
})


