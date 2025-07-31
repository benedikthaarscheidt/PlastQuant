library(testthat)
library(PlastQuant)

# Two-level simple case without covariates
test_that("simple two-level no covariates returns correct PILSM", {
  trait <- c(1, 3, 2, 4)
  env   <- c("A", "A", "B", "B")
  # LSMs equal means: A=2, B=3 → (3-2)/3 = 1/3
  expect_equal(calculate_PILSM(trait, env), 1/3)
})

# Three-level case without covariates
test_that("three-level no covariates returns correct PILSM", {
  trait <- c(1, 3, 2, 6, 4, 8)
  env   <- c("A", "A", "B", "B", "C", "C")
  # Means: A=2, B=4, C=6 → (6-2)/6 = 2/3
  expect_equal(calculate_PILSM(trait, env), 2/3)
})

# env = NULL: each obs its own level
test_that("env NULL treats each obs as unique level", {
  trait <- c(5, 10, 15)
  # LSMs equal trait: max=15, min=5 → (15-5)/15 = 2/3
  expect_equal(suppressWarnings(calculate_PILSM(trait, env = NULL)), 2/3)
})



# Zero maximum LSM warns and returns NA
test_that("zero maximum LSM returns NA with warning", {
  trait <- c(0, 0, 0, 0)
  env   <- c("A", "A", "B", "B")
  expect_warning(
    res <- calculate_PILSM(trait, env),
    "Maximum LSM is zero"
  )
  expect_true(is.na(res))
})

# Numeric vs. character env equivalence
test_that("numeric vs. character env gives same result", {
  trait <- c(1, 3, 2, 4)
  env_num <- c(1, 1, 2, 2)
  env_chr <- c("1", "1", "2", "2")
  expect_equal(
    calculate_PILSM(trait, env_num),
    calculate_PILSM(trait, env_chr)
  )
})

# Covariates length mismatch triggers error
test_that("covariate length mismatch triggers error", {
  trait <- c(1, 2, 3, 4)
  env   <- c("A", "A", "B", "B")
  cov   <- c(1, 2)  # length 2 != 4
  expect_error(
    calculate_PILSM(trait, env, covariates = cov),
    "must match length"
  )
})

# Non-vector/data.frame covariates triggers error
test_that("non-vector/data.frame covariates triggers error", {
  trait <- c(1, 2, 3)
  env   <- c("A", "B", "C")
  cov   <- matrix(1:3, ncol = 1)
  expect_error(
    calculate_PILSM(trait, env, covariates = cov),
    "covariates must be a vector or data.frame"
  )
})
