library(testthat)
library(PlastQuant)
test_that("errors on bad trait_values type or length", {
  expect_error(calculate_TPS("foo", 1:2, "A", "B"),
               "`trait_values` must be a numeric vector")
  expect_error(calculate_TPS(numeric(0), c(), "A", "B"),
               "`trait_values` must have length >= 1")
})

test_that("errors on missing or mismatched env_values", {
  tv <- 1:3
  expect_error(calculate_TPS(tv, NULL, "A", "B"),
               "`env_values` must be provided")
  expect_error(calculate_TPS(tv, 1:2, "A", "B"),
               "same length as `trait_values`")
})

test_that("errors on unknown native or transplanted envs", {
  tv  <- 1:4
  env <- rep(c("X","Y"), each = 2)
  expect_error(calculate_TPS(tv, env, "Z", "Y"),
               "`native_env` not found")
  expect_error(calculate_TPS(tv, env, "X", "Z"),
               "`transplanted_env` not found")
})

test_that("basic TPS calculation works", {
  tv  <- c(1, 2)
  env <- c("A", "B")
  res <- calculate_TPS(tv, env, "A", "B")
  expect_equal(res$raw, 1)
  expect_equal(res$overall, 1)
})

test_that("factor and character env_values work identically", {
  tv  <- c(2, 4, 3, 6)
  env <- rep(c("A","B"), each = 2)
  r1  <- calculate_TPS(tv, as.factor(env), "A", "B")
  r2  <- calculate_TPS(tv, as.character(env), "A", "B")
  expect_equal(r1$raw, r2$raw)
  expect_equal(r1$overall, r2$overall)
})

test_that("drops zeros when drop_zero = TRUE", {
  tv <- c(0, 5, 1, 4)
  # NOTE: use 'times = 2' so that (0,1) are A's and (5,4) are B's
  env <- rep(c("A","B"), times = 2)
  expect_warning(
    res <- calculate_TPS(tv, env, "A", "B", drop_zero = TRUE),
    "Dropping 1 zero native value"
  )
  expect_equal(res$raw, 3)      # (4-1)/1
  expect_equal(res$overall, 3)
})

test_that("Inf if drop_zero = FALSE, no warning", {
  tv  <- c(0, 5)
  env <- c("A", "B")
  expect_silent(
    res <- calculate_TPS(tv, env, "A", "B", drop_zero = FALSE)
  )
  expect_equal(res$raw, Inf)
  expect_equal(res$overall, Inf)
})

test_that("returns NA when all dropped", {
  tv  <- c(0, 0)
  env <- c("A", "B")
  expect_warning(r <- calculate_TPS(tv, env, "A", "B", drop_zero = TRUE),
                 "No valid pairs remain")
  expect_equal(length(r$raw), 0)
  expect_true(is.na(r$overall))
})

test_that("custom summary_fun and na.rm work", {
  tv  <- c(1, NA, 2, NA)
  env <- rep(c("A","B"), each = 2)
  # raw: (2-1)/1 = 1; (NA-NA)/NA = NaN
  r1 <- calculate_TPS(tv, env, "A","B", summary_fun = mean, na.rm = TRUE)
  expect_equal(r1$overall, 1)
  r2 <- calculate_TPS(tv, env, "A","B", summary_fun = median, na.rm = TRUE)
  expect_equal(r2$overall, 1)
})
