library(testthat)
library(PlastQuant)

test_that("non-numeric trait_values errors", {
  expect_error(calculate_RPI("a"), "`trait_values` must be a numeric vector")
})

test_that("too few observations returns NA warning", {
  expect_warning(res <- calculate_RPI(c(NA, 5)), "Not enough non-NA observations")
  expect_true(is.na(res))
})

test_that("env length mismatch errors", {
  expect_error(calculate_RPI(1:3, env_values = 1:2), "must have the same length")
})


test_that("less than two unique envs returns NA warning", {
  expect_warning(res <- calculate_RPI(1:3, env_values = c(1,1,1)), "Less than two unique environments")
  expect_true(is.na(res))
})

test_that("specific env pair computes correctly", {
  trait <- c(1,2,3,4)
  env   <- c("A","A","B","B")
  # env1="A"->trait c(1,2); env2="B"->trait c(3,4); min_len=2; ratios=c(2/5,2/6)=c(0.4,0.3333)
  exp <- mean(c(abs(3-1)/(3+1), abs(4-2)/(4+2)))
  expect_equal(calculate_RPI(trait, env, env1="A", env2="B"), exp)
})

test_that("all pairs computes grand mean correctly", {
  trait <- c(1,2,3)
  env   <- c(1,2,3)
  # pairs: (1,2)->1/3; (1,3)->2/4=0.5; (2,3)->1/5; mean=(1/3+0.5+0.2)/3
  exp <- mean(c(1/3,0.5,1/5))
  expect_equal(calculate_RPI(trait, env), exp)
})

test_that("division by zero yields NA in ratio", {
  trait <- c(1,-1)
  env   <- c(1,2)
  # denom=0 -> NA, mean(NA)=NA
  expect_true(is.na(calculate_RPI(trait, env)))
})
