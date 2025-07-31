library(testthat)
library(PlastQuant)

test_that("calculate_reaction_norm_slope returns NA for vectors of length < 2", {
  expect_true(is.na(calculate_reaction_norm_slope(numeric(0))))
  expect_true(is.na(calculate_reaction_norm_slope(42)))
})

test_that("calculate_reaction_norm_slope uses sequential environments by default", {
  trait <- c(5, 6, 7, 8, 9)
  expect_equal(calculate_reaction_norm_slope(trait), 1)
  trait2 <- c(2, 4, 6, 8)
  expect_equal(calculate_reaction_norm_slope(trait2), 2)
})

test_that("calculate_reaction_norm_slope accepts custom environments", {

  env <- c(10, 20, 30, 40)
  trait <- c(1, 3, 5, 7)
  expect_equal(calculate_reaction_norm_slope(trait, environments = env), 0.2)

  # non–unit spacing example
  env2 <- c(0, 5, 15)
  trait2 <- c(0, 10, 30)
  # slope = (30 - 0)/(15 - 0) = 2
  expect_equal(calculate_reaction_norm_slope(trait2, environments = env2), 2)
})

test_that("calculate_reaction_norm_slope handles non‐numeric or NA inputs via lm defaults", {
  # lm will error if there are NAs in the response
  expect_error(calculate_reaction_norm_slope(c(1, NA, 3)))

  # If env contains NA, lm errors
  expect_error(
    calculate_reaction_norm_slope(c(1, 2, 3), environments = c(1, NA, 3))
  )
})
