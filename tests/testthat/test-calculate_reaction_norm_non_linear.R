library(testthat)
library(PlastQuant)  # replace with your actual package name

test_that("errors on non-numeric trait_values", {
  expect_error(
    calculate_reaction_norm_non_linear(letters[1:5]),
    "`trait_values` must be a numeric vector"
  )
})

test_that("errors on NA in trait_values", {
  expect_error(
    calculate_reaction_norm_non_linear(c(1, NA, 3)),
    "contains missing values"
  )
})

test_that("errors on invalid degree", {
  expect_error(
    calculate_reaction_norm_non_linear(1:5, degree = 0),
    "`degree` must be a single integer"
  )
  expect_error(
    calculate_reaction_norm_non_linear(1:5, degree = 2.5),
    "`degree` must be a single integer"
  )
})

test_that("errors on bad environments", {
  expect_error(
    calculate_reaction_norm_non_linear(1:5, environments = letters[1:5]),
    "`environments` must be a numeric vector"
  )
  expect_error(
    calculate_reaction_norm_non_linear(1:5, environments = c(1, 2, NA, 4, 5)),
    "contains missing values"
  )
  expect_error(
    calculate_reaction_norm_non_linear(1:5, environments = 1:4),
    "same length as `trait_values`"
  )
})
