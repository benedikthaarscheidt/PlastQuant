library(testthat)
library(PlastQuant)

library(testthat)

context("calculate_FPI")

test_that("errors on bad trait_values and length < 2", {
  expect_error(calculate_FPI("foo"),
               "`trait_values` must be a numeric vector")
  expect_error(calculate_FPI(numeric(1)),
               "Need at least two trait values")
})

test_that("errors on mismatched env_values length", {
  tv <- rnorm(3)
  expect_error(calculate_FPI(tv, env_values = 1:2),
               "`env_values` must have same length")
})

test_that("default env_values gives correct FPI for two points", {
  res <- calculate_FPI(c(1, 2))
  expect_equal(res$pairwise, setNames(1, "1–2"))
  expect_equal(res$overall, 1)
})

test_that("numeric vs factor env_values give same values", {
  tv  <- c(1, 3, 5, 7)
  env <- rep(c("X", "Y"), each = 2)
  res_num <- calculate_FPI(tv, env_values = as.numeric(factor(env)))
  res_fac <- calculate_FPI(tv, env_values = factor(env))
  expect_equal(unname(res_num$pairwise),
               unname(res_fac$pairwise))
  expect_equal(res_num$overall, res_fac$overall)
})

test_that("character env_values also work", {
  tv   <- c(10, 15)
  res1 <- calculate_FPI(tv, env_values = c("A","B"))
  expect_equal(res1$overall, 0.5)
})

test_that("pairwise mean matches hand calculation", {
  tv  <- c(1,3,5,7)
  env <- rep(c("X","Y"), each = 2)
  # (5-1)/1 = 4; (7-3)/3 = 4/3; mean = (4 + 4/3)/2 = 8/3
  res <- calculate_FPI(tv, env)
  expect_equal(unname(res$pairwise), 8/3)
  expect_equal(res$overall, 8/3)
})

test_that("control_stress numeric 0/1 gives correct result", {
  tv <- c(2,4)
  cs <- c(0,1)
  res <- calculate_FPI(tv, control_stress = cs)
  expect_equal(unname(res$pairwise), 1)
  expect_equal(res$overall, 1)
})

test_that("control_stress logical is coerced", {
  tv <- c(2,8)
  cs <- c(FALSE, TRUE)
  res <- calculate_FPI(tv, control_stress = cs)
  expect_equal(res$overall, 3)
})

test_that("control_stress 'Control'/'Stress' works", {
  tv <- c(5,10)
  cs <- c("Control","Stress")
  res <- calculate_FPI(tv, control_stress = cs)
  expect_equal(res$overall, 1)
})

test_that("drop_zero = FALSE yields Inf, and no warning", {
  tv <- c(0, 5)
  cs <- c(0, 1)
  expect_silent(
    res <- calculate_FPI(tv, control_stress = cs, drop_zero = FALSE)
  )
  expect_equal(res$overall, Inf)
})



test_that("custom summary_fun and na.rm behaviour", {
  tv  <- c(1, NA, 3, 5)
  env <- rep(c("A","B"), each = 2)
  # pair (3-1)/1 = 2; (5-NA)/NA = NaN
  res1 <- calculate_FPI(tv, env, summary_fun = mean, na.rm = TRUE)
  expect_equal(res1$overall, 2)
  res2 <- calculate_FPI(tv, env, summary_fun = median, na.rm = TRUE)
  expect_equal(res2$overall, 2)
})

test_that("single unique environment warns and returns NA", {
  tv  <- c(1,2,3)
  env <- rep("Z", 3)
  expect_warning(
    res <- calculate_FPI(tv, env_values = env),
    "Only one unique environment"
  )
  expect_true(is.na(res$overall))
  expect_true(is.na(res$pairwise))
})
