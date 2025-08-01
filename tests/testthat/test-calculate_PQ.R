library(testthat)
library(PlastQuant)

test_that("basic PQ with equidistant env computes correctly", {
  trait <- c(1, 3, 6)
  # env=1,2,3 => trait_range=5, env_range=2 => PQ=2.5
  expect_equal(calculate_PQ(trait), 5/2)
})

test_that("provided numeric env_values works correctly", {
  trait <- c(2, 5, 9)
  env   <- c(10, 20, 40)
  # trait_range=7, env_range=30 => PQ=7/30
  expect_equal(calculate_PQ(trait, env), 7/30)
})

test_that("factor env_values coercion", {
  trait <- c(1,4,7)
  env_f  <- factor(c("1","3","6"))
  expect_equal(calculate_PQ(trait, env_f), (7-1)/(6-1))
})

test_that("character env_values non-numeric errors", {
  trait <- c(1,2,3)
  expect_error(calculate_PQ(trait, c("a","b","c")), "Cannot coerce character `env_values` to numeric.")
})

test_that("insufficient trait_values errors", {
  expect_error(calculate_PQ(1), "At least two non-NA trait values")
})

test_that("env length mismatch errors", {
  expect_error(calculate_PQ(c(1,2,3), c(1,2)), "same length as `trait_values`")
})

test_that("zero env range returns NA with warning", {
  trait <- c(1,2,3)
  env   <- c(5,5,5)
  expect_warning(res <- calculate_PQ(trait, env), "Environmental range is zero")
  expect_true(is.na(res))
})
