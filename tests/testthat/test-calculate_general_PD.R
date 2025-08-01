library(testthat)
library(PlastQuant)

test_that("trait_values input validation", {
  expect_error(calculate_general_PD("a", env_values = 1:2), "'trait_values' must be numeric")
  expect_error(calculate_general_PD(1, env_values = 1), "At least two trait values are required")
})

test_that("env_values default and length check", {
  # default env_values should produce no error
  expect_silent(calculate_general_PD(c(1,2)))
  # mismatched lengths
  expect_error(calculate_general_PD(1:3, env_values = 1:2), "'env_values' must be the same length as 'trait_values'.")
})

test_that("NA filtering returns NA with warning when insufficient data", {
  expect_warning(res <- calculate_general_PD(c(NA, NA, 3)), "Not enough valid data points; returning NA")
  expect_true(is.na(res))
})

test_that("variability method returns max-min", {
  x <- c(5, 3, NA, 9)
  expect_equal(calculate_general_PD(x, env_values = 1:length(x), method = "variability"), 9 - 3)
})

test_that("reference method computes correct PD", {
  trait <- c(2,4,6,8)
  env   <- c("A","A","B","B")
  # means A=3, B=7; abs diffs from ref(A): B-A = 4
  expect_equal(
    calculate_general_PD(trait, env_values = env, method = "reference"),
    4
  )
  # with three envs
  trait2 <- c(1,3,5,7,9,11)
  env2   <- c(1,1,2,2,3,3)
  # means:1->2,2->6,3->10; ref=1 -> diffs 4 and 8 -> avg=6
  expect_equal(
    calculate_general_PD(trait2, env_values = env2, method = "reference"),
    6
  )
})

test_that("pairwise method computes correct average pairwise PD", {
  trait <- c(1,2,3,4)
  env   <- c(1,1,2,2)
  # pairwise envs: (1vs2): mean(abs((1,2)-(3,4))) = mean(c(2,2)) = 2
  expect_equal(
    calculate_general_PD(trait, env_values = env, method = "pairwise"),
    2
  )
  # with three envs
  trait3 <- c(1,2,3,4,5,6)
  env3   <- c("X","X","Y","Y","Z","Z")
  # means per pair: X-Y: mean(abs((1,2)-(3,4)))=mean(c(2,2))=2
  # X-Z: mean(abs((1,2)-(5,6)))=mean(c(4,4))=4
  # Y-Z: mean(abs((3,4)-(5,6)))=mean(c(2,2))=2
  # avg=(2+4+2)/3 = 8/3
  expect_equal(
    calculate_general_PD(trait3, env_values = env3, method = "pairwise"),
    8/3
  )
})

test_that("control_stress_vector method", {
  trait <- c(1,3,5,7)
  csv1  <- c(0,1,0,1)
  # pairs: (1,3)->2, (5,7)->2 -> mean=2
  expect_equal(
    calculate_general_PD(trait, control_stress_vector = csv1),
    2
  )
  csv2 <- c("Control","Stress","Control","Stress")
  expect_equal(
    calculate_general_PD(trait, control_stress_vector = csv2),
    2
  )
  # unequal lengths should error
  expect_error(
    calculate_general_PD(trait, control_stress_vector = c(0,1,0)),
    "must be the same length"
  )
})

test_that("invalid method errors", {
  expect_error(
    calculate_general_PD(1:4, env_values = 1:4, method = "foo"),
    "must be 'pairwise', 'reference', or 'variability'"
  )
})
