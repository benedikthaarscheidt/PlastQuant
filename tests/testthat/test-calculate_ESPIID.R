library(testthat)
library(PlastQuant)

# Non-numeric trait_values triggers error
test_that("non-numeric trait_values triggers error", {
  expect_error(
    calculate_espiid("a"),
    "`trait_values` must be a numeric vector."
  )
})

# use_median and aggregate must be single logical values
test_that("use_median must be logical scalar", {
  expect_error(
    calculate_espiid(1:3, use_median = "TRUE"),
    "`use_median` must be a single TRUE or FALSE value."
  )
})

test_that("aggregate must be logical scalar", {
  expect_error(
    calculate_espiid(1:3, aggregate = NA),
    "`aggregate` must be a single TRUE or FALSE value."
  )
})

# Too few non-NA observations yields warning and returns NA
test_that("too few non-NA observations returns NA with warning", {
  vec <- c(NA, NA, 5)
  expect_warning(
    res <- calculate_espiid(vec),
    "At least two non-NA trait values are required; ESPIID cannot be calculated."
  )
  expect_true(is.na(res))
})

# Basic aggregated ESPIID computation
test_that("basic aggregated ESPIID computes correctly", {
  # trait 1,3,6: differences: (1,3)=2/1=2; (1,6)=5/2=2.5; (3,6)=3/1=3; mean = (2+2.5+3)/3
  trait <- c(1,3,6)
  expected <- mean(c(2/1,5/2,3/1))
  expect_equal(
    calculate_espiid(trait),
    expected
  )
})

# aggregate=FALSE returns named vector of pairwise ESPIID
test_that("aggregate=FALSE returns named vector of ESPIID", {
  trait <- c(2, 5, 8)
  # pairs: 1-2:|5-2|/1=3; 1-3:|8-2|/2=3; 2-3:|8-5|/1=3
  res_vec <- calculate_espiid(trait, aggregate = FALSE)
  expect_named(res_vec, c("1-2", "1-3", "2-3"))
  expect_equal(unname(res_vec), c(3,3,3))
})

# use_median has no effect for single observations
test_that("use_median yields same as mean when aggregate=FALSE for single obs per env", {
  trait <- c(4, 7)
  vec_mean   <- calculate_espiid(trait, use_median = FALSE, aggregate = FALSE)
  vec_median <- calculate_espiid(trait, use_median = TRUE,  aggregate = FALSE)
  expect_equal(vec_mean, vec_median)
})

# Handling NAs in trait_values
test_that("handles NAs in trait_values correctly when aggregate=FALSE", {
  trait <- c(NA, 3, 7)
  # after dropping NA: trait=c(3,7), only pair 1-2 => |7-3|/1=4
  res <- calculate_espiid(trait, aggregate = FALSE)
  expect_equal(unname(res), 4)
})
