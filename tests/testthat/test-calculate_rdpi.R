library(testthat)
library(PlastQuant)

# Non-numeric trait_values triggers error
test_that("non-numeric trait_values triggers error", {
  expect_error(
    calculate_rdpi("a", env_values = c(1, 2)),
    "trait_values must be a numeric vector."
  )
})

# trait_values length < 2 triggers error
test_that("too few trait_values triggers error", {
  expect_error(
    calculate_rdpi(c(1)),
    "trait_values must contain at least two observations."
  )
})

# Mismatched lengths triggers error
test_that("mismatched lengths triggers error", {
  trait <- c(1, 2, 3)
  env   <- c("A", "B")
  expect_error(
    calculate_rdpi(trait, env_values = env),
    "trait_values and env_values must have the same length."
  )
})

# Only one environment level triggers error
test_that("single environment level triggers error", {
  trait <- c(1, 2, 3)
  env   <- c("A", "A", "A")
  expect_error(
    calculate_rdpi(trait, env_values = env),
    "At least two environment levels are required to calculate RDPI."
  )
})

# Basic two-environment case without env_values
test_that("basic RDPI with default env_values computes correctly", {
  tv <- c(1, 3)
  # default env_values=> unique levels => each unique => two levels, means=1,3 => RDPI=|3-1|/(3+1)=0.5
  expect_equal(
    calculate_rdpi(tv),
    0.5
  )
})

# Grouped environments compute based on means
test_that("grouped environments compute RDPI correctly", {
  tv  <- c(1, 3, 2, 4)
  env <- c("A", "A", "B", "B")
  # means A=2, B=3 => RDPI = (3-2)/(3+2) = 1/5
  expect_equal(
    calculate_rdpi(tv, env_values = env),
    1/5
  )
})

# Handling NA values within groups
test_that("handles NA values within groups correctly", {
  tv  <- c(NA, 3, 2, NA)
  env <- c("A", "A", "B", "B")
  # means A=3, B=2 => RDPI = |2-3|/(2+3)=1/5
  expect_equal(
    calculate_rdpi(tv, env_values = env),
    1/5
  )
})

# All-group NA leads to warning and NA result
test_that("all-group NA triggers warning and returns NA", {
  tv  <- c(NA_real_, NA_real_, 1, 2)
  env <- c("A", "A", "B", "B")
  # A group all NA => drop A => only one mean => warning + NA
  expect_warning(
    res <- calculate_rdpi(tv, env_values = env),
    "Not enough valid environment means"
  )
  expect_true(is.na(res))
})

# Zero denominators produce RDPI = 0
test_that("zero means produce RDPI zero without warning", {
  tv  <- c(0, 0, 0, 0)
  env <- c("A", "A", "B", "B")
  # means 0,0 => sums zero => each pair gives 0 => RDPI=0
  expect_equal(
    calculate_rdpi(tv, env_values = env),
    0
  )
})
