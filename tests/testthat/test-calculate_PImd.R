library(testthat)
library(PlastQuant)

# Simple two-level case: character env
test_that("simple two-level returns correct PImd", {
  trait <- c(1, 3, 2, 4)
  env   <- c("A", "A", "B", "B")
  # medians: A=2, B=3 → (3-2)/3 = 1/3
  expect_equal(calculate_PImd(trait, env), 1/3)
})

# Three-level case: factor env
test_that("three-level returns correct PImd", {
  trait <- c(1, 3, 2, 6, 4, 8)
  env   <- factor(c("A", "A", "B", "B", "C", "C"))
  # medians: A=2, B=4, C=6 → (6-2)/6 = 2/3
  expect_equal(calculate_PImd(trait, env), 2/3)
})

# env = NULL: each observation its own level
test_that("env NULL treats each obs as unique level", {
  trait <- c(5, 10, 15)
  # medians = same as values → max=15, min=5 → (15-5)/15 = 2/3
  expect_equal(calculate_PImd(trait, env = NULL), 2/3)
})

# NA median handling: some NAs but at least one median
test_that("handles partially NA trait_values correctly", {
  trait <- c(NA, 2, 4, NA)
  env   <- c("A", "A", "B", "B")
  # medians: A=2, B=4 → (4-2)/4 = 1/2
  expect_equal(calculate_PImd(trait, env), 1/2)
})

# Zero maximum median: warn and return NA
test_that("zero maximum median returns NA with warning", {
  trait <- c(0, 0, 0, 0)
  env   <- c("A", "A", "B", "B")
  expect_warning(res <- calculate_PImd(trait, env), "Maximum median is zero")
  expect_true(is.na(res))
})

# Length mismatch between trait_values and env
test_that("length mismatch triggers error", {
  trait <- c(1, 2, 3)
  env   <- c("A", "B")
  expect_error(
    calculate_PImd(trait, env),
    "must match length"
  )
})

# Non-numeric trait_values triggers error
test_that("non-numeric trait_values triggers error", {
  trait <- as.character(c(1, 2, 3))
  env   <- c("A", "A", "B")
  expect_error(
    calculate_PImd(trait, env),
    "trait_values must be a numeric vector"
  )
})

# Numeric env vs character env equivalence
test_that("numeric vs character env gives same result", {
  trait <- c(1, 3, 2, 4)
  env_num <- c(1, 1, 2, 2)
  env_chr <- c("1", "1", "2", "2")
  expect_equal(
    calculate_PImd(trait, env_num),
    calculate_PImd(trait, env_chr)
  )
})
