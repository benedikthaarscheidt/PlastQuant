library(testthat)
library(PlastQuant)

test_that("non-numeric trait_values errors", {
  expect_error(
    calculate_NRW("a"),
    "'trait_values' must be a numeric vector"
  )
})

test_that("across must be logical scalar", {
  expect_error(
    calculate_NRW(1:3, across = NA),
    "'across' must be a single TRUE or FALSE"
  )
})

test_that("env_values length mismatch errors", {
  expect_error(
    calculate_NRW(1:3, env_values = 1:2),
    "must be the same length"
  )
})


test_that("group_values length mismatch errors", {
  expect_error(
    calculate_NRW(1:3, group_values = 1:2),
    "'group_values' must be the same length as 'trait_values'"
  )
})

test_that("all NA values returns NA with warning", {
  expect_warning(
    res <- calculate_NRW(c(NA_real_, NA_real_)),
    "No valid data points after removing NAs"
  )
  expect_true(is.na(res))
})

test_that("basic NRW across computes correct range", {
  expect_equal(
    calculate_NRW(c(2, 5, 9), across = TRUE),
    9 - 2
  )
})

test_that("NRW across ignores env_values", {
  expect_equal(
    calculate_NRW(c(2, 6, 3), env_values = c(1, 2, 3), across = TRUE),
    6 - 2
  )
})

test_that("NRW per-environment computes correctly", {
  trait <- c(1, 3, 2, 5, 4)
  env   <- c("A", "A", "B", "B", "B")
  # mean_A = 2, mean_B = 11/3 => range = 11/3 - 2
  expected <- (11/3) - 2
  expect_equal(
    calculate_NRW(trait, env, across = FALSE),
    expected
  )
})

test_that("default env within returns zero for constant trait", {
  trait <- c(7, 7, 7)
  expect_equal(
    calculate_NRW(trait, across = FALSE),
    0
  )
})

test_that("NRW per group computes correctly", {
  trait <- c(1, 4, 2, 5, 3)
  env   <- c("A", "B", "A", "B", "B")
  grp   <- c("G1", "G1", "G2", "G2", "G2")
  # For G1: means A=1, B=4 => 4-1 = 3
  # For G2: means A=2, B=(5+3)/2 = 4 => 4-2 = 2
  res <- calculate_NRW(trait, env, group_values = grp, across = FALSE)
  expect_named(res, c("G1", "G2"))
  expect_equal(unname(res["G1"]), 3)
  expect_equal(unname(res["G2"]), 2)
})

test_that("character env_values coercion works", {
  trait <- c(2, 4, 6)
  env   <- c("X", "Y", "X")
  # mean_X = 4, mean_Y = 4 => range = 0
  expect_equal(
    calculate_NRW(trait, env, across = FALSE),
    0
  )
})
