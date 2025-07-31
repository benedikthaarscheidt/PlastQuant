library(testthat)
library(PlastQuant)

# Single-vector mode
test_that("vector input computes (max - min)/max", {
  v <- c(2, 4, 6, 8)
  expect_equal(calculate_Phenotypic_Plasticity_Index(v), 0.75)
})

# Vector with NAs
test_that("vector with NAs skips NA values", {
  v <- c(5, NA, 15, NA)
  expect_equal(calculate_Phenotypic_Plasticity_Index(v), 2/3)
})

# Single-column data frame returns scalar
test_that("data.frame with one named trait column returns scalar", {
  df <- data.frame(tr1 = c(1, 3, 5, 7))
  expect_equal(calculate_Phenotypic_Plasticity_Index(df, trait_cols = "tr1"), 6/7)
})

# Multiple trait columns return named vector
test_that("data.frame with multiple trait columns returns vector of Pis", {
  df <- data.frame(
    t1 = c(10, 20, 30),
    t2 = c(2,   5,  8)
  )
  pis <- calculate_Phenotypic_Plasticity_Index(df, trait_cols = c("t1", "t2"))
  expect_named(pis, c("t1", "t2"))
  expect_equal(unname(pis["t1"]), 2/3)
  expect_equal(unname(pis["t2"]), 0.75)
})

# Numeric indexing of columns works and returns scalar
test_that("numeric indexing of columns works and returns scalar", {
  df <- data.frame(a = c(1, 4), b = c(2, 8))
  expect_equal(calculate_Phenotypic_Plasticity_Index(df, trait_cols = 2), (8 - 2) / 8)
})

# Error on nonexistent character column
test_that("nonexistent character trait_cols triggers error", {
  df <- data.frame(x = 1:3)
  expect_error(
    calculate_Phenotypic_Plasticity_Index(df, trait_cols = "y"),
    "Columns not found: y"
  )
})

# Error on out-of-range numeric index
test_that("out-of-range numeric trait_cols triggers error", {
  df <- data.frame(x = 1:3)
  expect_error(
    calculate_Phenotypic_Plasticity_Index(df, trait_cols = 5),
    "Numeric trait_cols indices out of range"
  )
})
