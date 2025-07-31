library(testthat)
library(PlastQuant)

# Basic two-level comparison without covariates
test_that("simple two-level no covariate returns 100", {
  trait <- c(10, 10, 20, 20)
  env   <- c("A", "A", "B", "B")
  capture.output(ppf <- calculate_PPF(trait, env))
  expect_equal(unname(ppf), 100)
})

# Covariates perfectly confounded must trigger a warning and return unadjusted PPF
test_that("single covariate vector perfectly confounded triggers warning and returns 100", {
  trait <- c(10, 10, 20, 20)
  env   <- c("A", "A", "B", "B")
  cov   <- c(1, 1, 2, 2)
  expect_warning(
    capture.output(ppf_cov <- calculate_PPF(trait, env, cov)),
    "Collinearity detected"
  )
  expect_equal(unname(ppf_cov), 100)
})

# Multiple covariates perfectly confounded
test_that("multiple covariate data.frame perfectly confounded triggers warning and returns 100", {
  trait <- c(10, 10, 20, 20)
  env   <- c("A", "A", "B", "B")
  cov_df <- data.frame(x = c(1, 1, 2, 2), y = c(1, 1, 2, 2))
  expect_warning(
    capture.output(ppf_df <- calculate_PPF(trait, env, cov_df)),
    "Collinearity detected"
  )
  expect_equal(unname(ppf_df), 100)
})

# Three environments: check averaging over all environment pairs
test_that("three environments averages over all pairs", {
  trait <- c(10, 10, 20, 20, 30, 30)
  env   <- c("A", "A", "B", "B", "C", "C")
  expected <- mean(c(100, 200, 50))
  capture.output(ppf_all <- calculate_PPF(trait, env))
  expect_equal(unname(ppf_all), expected)
})

# Specific env_pairs argument
test_that("specific env_pairs argument works for one pair", {
  trait <- c(10, 10, 20, 20, 30, 30)
  env   <- c("A", "A", "B", "B", "C", "C")
  capture.output(ppf_AC <- calculate_PPF(trait, env, env_pairs = list(c("A", "C"))))
  expect_equal(unname(ppf_AC), 200)
})

# Error for single-level input
test_that("single-level environments throws an error", {
  expect_error(calculate_PPF(1:5, rep("X", 5)))
})
