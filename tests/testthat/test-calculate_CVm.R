library(testthat)
library(PlastQuant)

test_that("list‐of‐groups mode computes correctly", {
  trait_groups <- list(
    Env1 = c(10, 12, 11),
    Env2 = c(20, 22, 21),
    Env3 = c(16, 18, 17)
  )
  # group means = c(11, 21, 16); sd = sd(c(11,21,16)); mean = mean(c(11,21,16))
  expected <- sd(c(11, 21, 17)) / mean(c(11, 21, 17))
  expect_equal(calculate_CVm(trait_values = trait_groups), expected)
})

test_that("vector+labels mode matches list‐mode", {
  traits <- c(10,12,11, 20,22,21, 15,18,17)
  labels <- rep(c("Env1","Env2","Env3"), each = 3)
  expect_equal(
    calculate_CVm(trait_values = traits, group_labels = labels),
    calculate_CVm(trait_values = list(Env1=c(10,12,11), Env2=c(20,22,21), Env3=c(15,18,17)))
  )
})

test_that("errors if group_labels missing or mismatched in vector mode", {
  expect_error(
    calculate_CVm(trait_values = 1:5),
    "`group_labels` must be provided"
  )
  expect_error(
    calculate_CVm(trait_values = 1:4, group_labels = c("A","B","C")),
    "must have the same length"
  )
})

test_that("handles NA within groups (na.rm = TRUE)", {
  # Env1 mean = mean(c(1,NA,3), na.rm=TRUE)=2; Env2 mean = 4
  groups <- list(A = c(1, NA, 3), B = c(4, NA))
  expected_means <- c(A = mean(c(1,NA,3), na.rm = TRUE), B = mean(c(4,NA), na.rm = TRUE))
  expected <- sd(expected_means) / mean(expected_means)
  expect_equal(calculate_CVm(groups), expected)
})

test_that("zero between‐group variance returns 0", {
  # All group means identical → sd_means = 0 → CVm = 0
  groups <- list(G1 = c(5,5), G2 = c(5,5), G3 = c(5,5))
  expect_equal(calculate_CVm(groups), 0)
})

test_that("zero grand mean returns NA", {
  # Group means sum to zero → mean_means = 0 → returns NA to avoid div by zero
  groups <- list(G1 = c(-1,1), G2 = c(-2,2))
  # group_means = c(0, 0) → mean_means = 0
  expect_true(is.na(calculate_CVm(groups)))
})
