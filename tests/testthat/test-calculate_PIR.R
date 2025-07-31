library(testthat)
library(PlastQuant)

# Non-numeric trait_values triggers error
test_that("non-numeric trait_values triggers error", {
  expect_error(
    calculate_PIR("a", env_values = c(1,2,3)),
    "trait_values must be a numeric vector"
  )
})

# env_values length mismatch triggers error
test_that("env_values length mismatch triggers error", {
  expect_error(
    calculate_PIR(trait_values = 1:3, env_values = 1:2),
    "env_values must have the same length"
  )
})

# All observations missing returns NA with warning
test_that("all missing returns NA with warning", {
  tv <- c(NA_real_, NA_real_, NA_real_)
  ev <- c("A","B","C")  # still numeric check means env factor conversion
  expect_warning(
    res <- calculate_PIR(tv, ev),
    "All observations are missing; PIR cannot be calculated"
  )
  expect_true(is.na(res))
})

# Less than two environments triggers error
test_that("single environment triggers error", {
  tv <- c(1,2,3)
  ev <- c("A","A","A")
  expect_error(
    calculate_PIR(tv, ev),
    "At least two environments are required to calculate PIR"
  )
})

# Basic calculation without rgr_values
test_that("basic PIR without rgr_values computes correctly", {
  # Two environments: means A=2, B=4, rgr computed: NA, (4-2)/2=1; max rgr index=2
  tv <- c(1,3,3,5)
  ev <- c("A","A","B","B")
  # means: A=(1+3)/2=2, B=(3+5)/2=4; rgr: NA, (4-2)/2=1; idx=2, mean_at_max=4
  # PIR = (max_mean - min_mean)/mean_at_max = (4-2)/4 = 0.5
  expect_equal(unname(calculate_PIR(tv, ev)), 0.5)
})

# Calculation with provided rgr_values
test_that("PIR with provided rgr_values computes correctly", {
  tv <- c(2,4,6,8)
  ev <- c("A","A","B","B")
  # user-supplied RGR per observation; aggregated: A=mean(c(0.1,0.2))=0.15, B=mean(c(0.3,0.4))=0.35
  rgr <- c(0.1,0.2,0.3,0.4)
  # means: A=3, B=7; max_rgr idx=2 -> mean_at_max=7; PIR=(7-3)/7=4/7
  expect_equal(unname(calculate_PIR(tv, ev, rgr_values = rgr)), (7-3)/7)
})

# Default RGR all-denominator-zero yields NA with warning
test_that("default RGR all NA returns NA with warning", {
  tv <- c(0,0,5,5)
  ev <- c("A","A","B","B")
  expect_warning(
    res <- calculate_PIR(tv, ev),
    "All RGR values are NA; PIR cannot be calculated"
  )
  expect_true(is.na(res))
})


