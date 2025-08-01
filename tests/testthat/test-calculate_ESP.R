library(testthat)
library(PlastQuant)


test_that("non-numeric trait_values errors", {
  expect_error(
    calculate_ESP("a"),
    "'trait_values' must be a numeric vector"
  )
})



test_that("env_values length mismatch errors", {
  expect_error(
    calculate_ESP(1:3, env_values = 1:2),
    "must be the same length"
  )
})

test_that("env_subset wrong type errors", {
  expect_error(
    calculate_ESP(1:3, env_values = 1:3, env_subset = list(1,2)),
    "'env_subset' must be numeric, factor, or character"
  )
})

test_that("env_subset with no matches returns NA with warning", {
  expect_warning(
    res <- calculate_ESP(1:3, env_values = c("A","B","C"), env_subset = "Z"),
    "None of the 'env_subset' values match"
  )
  expect_true(is.na(res))
})

test_that("all NA trait_values returns NA with warning", {
  expect_error(
    res <- calculate_ESP(c(NA, NA, NA)),
    "'trait_values' must be a numeric vector."
  )
})

test_that("basic ESP calculation is correct", {
  trait <- c(1, 3, 5, 7)
  env   <- c(1, 1, 2, 2)
  # mean_all = 4; mean_env1 = 2 => -0.5; mean_env2 = 6 => +0.5; sum(abs) = 1.0
  expect_equal(
    calculate_ESP(trait, env_values = env),
    1.0
  )
})

test_that("env_subset filters environments correctly", {
  trait <- c(1, 3, 5, 7)
  env   <- c("X", "X", "Y", "Y")
  # full: 1.0
  expect_equal(
    calculate_ESP(trait, env_values = env),
    1.0
  )
  # subset to Y only: (6-4)/4 = 0.5
  expect_equal(
    calculate_ESP(trait, env_values = env, env_subset = "Y"),
    0.5
  )
})

test_that("constant trait across env returns zero", {
  trait <- rep(5, 6)
  env   <- rep(1:3, each = 2)
  expect_equal(
    calculate_ESP(trait, env_values = env),
    0
  )
})

test_that("single NA in an environment is ignored without warning", {
  trait <- c(NA, 2, 4, 6)
  env   <- c("A",  "A", "B", "B")
  # After removing the NA, env A still has one value (2) and env B has (4,6):
  # mean_all = (2+4+6)/3 = 4
  # ESP_A = (2-4)/4 = -0.5; ESP_B = (5-4)/4 = 0.25; sum(abs) = 0.75
  expect_silent({
    res <- calculate_ESP(trait, env_values = env)
  })
  expect_equal(res, 0.75)
})
