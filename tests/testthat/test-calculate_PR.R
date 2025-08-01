library(testthat)
library(PlastQuant)

test_that("non-numeric trait_values errors", {
  expect_error(calculate_PR("a"), "`trait_values` must be a numeric vector")
})

test_that("across must be logical scalar", {
  expect_error(calculate_PR(1:3, across = NA), "`across` must be a single TRUE or FALSE value")
})

test_that("env length mismatch errors", {
  expect_error(calculate_PR(1:3, env_values = 1:2), "must have the same length as `trait_values`")
})

test_that("all NA values returns NA with warning", {
  expect_warning(res <- calculate_PR(c(NA_real_, NA_real_)), "Not enough non-NA observations")
  expect_true(is.na(res))
})

test_that("basic PR across computes correct range", {
  expect_equal(calculate_PR(c(2,5,9)), 9-2)
})

test_that("PR across ignores env_values", {
  expect_equal(calculate_PR(c(2,6,3), env_values = c(1,2,3)), 6-2)
})

test_that("PR within environments computes correctly", {
  trait <- c(1,3,2,5,4)
  env   <- c("A","A","B","B","B")
  res <- calculate_PR(trait, env, across = FALSE)
  expect_named(res, c("A","B"))
  expect_equal(unname(res["A"]), 3-1)
  expect_equal(unname(res["B"]), 5-2)
})

test_that("default env within returns zeros for each obs", {
  trait <- c(7,7,7)
  res <- calculate_PR(trait, across = FALSE)
  expect_equal(unname(res), c(0,0,0))
})

test_that("character env_values coercion works", {
  trait <- c(2,4,6)
  env   <- c("X","Y","X")
  res <- calculate_PR(trait, env, across = FALSE)
  expect_named(res, c("X","Y"))
  expect_equal(unname(res["X"]), 6-2)
  expect_equal(unname(res["Y"]), 4-4)
})

