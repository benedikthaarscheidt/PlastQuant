library(testthat)
library(PlastQuant)  # replace with your actual package name

test_that("list‐of‐groups mode computes correctly", {
  trait_groups <- list(
    Env1 = c(10, 12, 11),
    Env2 = c(20, 22, 21),
    Env3 = c(16, 18, 17)
  )
  # group medians = c(11, 21, 17); sd = sd(c(11,21,17)); mean = mean(c(11,21,17))
  expected <- sd(c(11, 21, 17)) / mean(c(11, 21, 17))
  expect_equal(calculate_CVmd(trait_values = trait_groups), expected)
})

test_that("vector+labels mode matches list‐mode", {
  traits <- c(10,12,11, 20,22,21, 15,18,17)
  labels <- rep(c("Env1","Env2","Env3"), each = 3)
  expect_equal(
    calculate_CVmd(trait_values = traits, group_labels = labels),
    calculate_CVmd(trait_values = list(Env1=c(10,12,11),
                                       Env2=c(20,22,21),
                                       Env3=c(15,18,17)))
  )
})

test_that("errors if group_labels missing or mismatched in vector mode", {
  expect_error(
    calculate_CVmd(trait_values = 1:5),
    "`group_labels` must be provided"
  )
  expect_error(
    calculate_CVmd(trait_values = 1:4, group_labels = c("A","B","C")),
    "must have the same length"
  )
})

test_that("handles NA within groups (na.rm = TRUE)", {
  # Env1 median = median(c(1,NA,3), na.rm=TRUE)=2; Env2 median = median(c(4,NA), na.rm=TRUE)=4
  groups <- list(A = c(1, NA, 3), B = c(4, NA))
  medians <- c(A = median(c(1,NA,3), na.rm = TRUE), B = median(c(4,NA), na.rm = TRUE))
  expected <- sd(medians) / mean(medians)
  expect_equal(calculate_CVmd(groups), expected)
})

test_that("zero between‐group variance returns 0", {
  # All group medians identical → sd_med = 0 → CVmd = 0
  groups <- list(G1 = c(5,5), G2 = c(5,5), G3 = c(5,5))
  expect_equal(calculate_CVmd(groups), 0)
})

test_that("zero grand mean returns NaN", {
  # Group medians sum to zero → mean_med = 0 → result = 0/0 = NaN
  groups <- list(G1 = c(-1,1), G2 = c(-2,2))
  out <- calculate_CVmd(groups)
  expect_true(is.na(out))
})
