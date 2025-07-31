library(testthat)
library(PlastQuant)

test_that("returns NA for fewer than 2 values", {
  expect_true(is.na(calculate_RC(numeric(0))))
  expect_true(is.na(calculate_RC(42)))
})

test_that("errors on missing values in input", {
  expect_error(
    calculate_RC(c(1, NA, 3)),
    "missing values"
  )
})

test_that("errors when fractions out of [0,1]", {
  expect_error(
    calculate_RC(1:5, lower_fraction = -0.1, upper_fraction = 0.5),
    "must be between 0 and 1"
  )
  expect_error(
    calculate_RC(1:5, lower_fraction = 0.5, upper_fraction = 1.1),
    "must be between 0 and 1"
  )
})

test_that("errors when lower_fraction >= upper_fraction", {
  expect_error(
    calculate_RC(1:5, lower_fraction = 0.8, upper_fraction = 0.2),
    "must be less than"
  )
})

test_that("correct RC on simple sorted data", {
  # trait = 1:10, length=10
  # default lower=0.5 → floor(10*0.5)=5 → lower=1:5, mean_low=3
  # default upper=0.5 → upper_boundary_index = ceiling(10*(1-0.5))=ceiling(5)=5 → upper=(6:10), mean_high=8
  # RC = 8/3
  expect_equal(
    calculate_RC(1:10, lower_fraction = 0.5, upper_fraction = 0.5),
    8/3
  )
})

test_that("correct RC with unsorted input", {
  vec <- seq(1:10)
  # sorted = 10,20,30,50,80,100 (n=6)
  # lower_fraction=0.3 → floor(6*0.3)=1 → lower=10, mean_low=10
  # upper_fraction=0.2 → upper_boundary_index=ceiling(6*0.8)=ceiling(4.8)=5 → upper=(6th index 100), mean_high=100
  expect_equal(
    calculate_RC(vec, lower_fraction = 0.3, upper_fraction = 0.8),
    9.5/2
  )
})

test_that("returns NA when mean_low or mean_high is zero", {
  # mean_low zero
  vec1 <- c(0, 0, 5, 10)
  # lower_fraction=0.5 → floor(4*0.5)=2 → lower=0,0 → mean_low=0 → NA
  expect_true(is.na(calculate_RC(vec1, lower_fraction = 0.5, upper_fraction = 0.75)))

  # mean_high zero
  vec2 <- c(1, 2, 0, 0)
  # sorted=0,0,1,2, lower=0,0 → mean_low=0 (NA triggered there already)
  expect_true(is.na(calculate_RC(vec2, lower_fraction = 0.5, upper_fraction = 0.75)))
})


