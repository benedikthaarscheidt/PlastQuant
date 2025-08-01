library(testthat)
library(PlastQuant)

# Basic two-point case: equidistant envs
test_that("simple ESPI with equidistant env computes correctly", {
  trait <- c(1, 3)
  # env = c(1,2): means are trait themselves; max_mean=3, min_mean=1; range=
  expect_equal(calculate_ESPI(trait), (3-1)/(1))
})

# Provided env_values numeric
test_that("ESPI with provided numeric env_values works", {
  trait <- c(2, 4, 6)
  env   <- c(10, 20, 40)
  # means = trait; max_mean=6,min_mean=2; env_range=40-10=30
  expect_equal(calculate_ESPI(trait, env), (6-2)/30)
})

# Factor env_values coercion
test_that("ESPI coerces factor env_values correctly", {
  trait <- c(5, 7, 9)
  env_fac <- factor(c("1","2","4"))
  expect_equal(calculate_ESPI(trait, env_fac), (9-5)/(3))
})

# Character env_values coercion and error on non-numeric
test_that("character env_values numeric coercion and error", {
  trait <- c(1,2,3)
  env_chr <- c("1","2","a")
  expect_error(calculate_ESPI(trait, env_chr), "Cannot coerce character `env_values` to numeric.")
})

# Missing observations dropped
test_that("ESPI drops NA observations and computes", {
  trait <- c(NA, 4, 6)
  env   <- c(1, 2, 3)
  # drop first: trait=c(4,6), env=c(2,3), max_mean=6,min_mean=4,range=1
  expect_equal(calculate_ESPI(trait, env), (6-4)/(3-2))
})

# All observations NA triggers warning and NA result
test_that("all NA trait_values returns NA with warning", {
  trait <- c(NA_real_, NA_real_)
  env   <- c(1,2)
  expect_warning(res <- calculate_ESPI(trait, env), "All observations are NA")
  expect_true(is.na(res))
})

# Zero environmental range triggers warning and NA result
test_that("zero env range returns NA with warning", {
  trait <- c(1,2,3)
  env   <- c(5,5,5)
  expect_warning(res <- calculate_ESPI(trait, env), "Environmental range is zero")
  expect_true(is.na(res))
})

# Input validation errors: trait non-numeric, length mismatch, too few obs
test_that("error on non-numeric trait_values", {
  expect_error(calculate_ESPI("a", 1:3), "`trait_values` must be a numeric vector")
})

test_that("error on trait length < 2", {
  expect_error(calculate_ESPI(1), "must contain at least two observations")
})

test_that("error on env length mismatch", {
  expect_error(calculate_ESPI(c(1,2), c(1)), "`env_values` must have the same length")
})
