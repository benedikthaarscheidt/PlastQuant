
library(testthat)
library(PlastQuant)



df_test4 <- data.frame(
  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),
  Column1 = c(rep(2, 10), rep(3, 10), rep(1, 10)),
  Column2 = c(rep(1, 10), rep(1, 10), rep(1, 10))
)
df_test5 <- data.frame(
  Column0 = c(rep(3, 10), rep(2, 10), rep(1, 10)),
  Column1 = c(rep(3, 10), rep(2, 10), rep(1, 10)),
  Column2 = c(rep(3, 10), rep(2, 10), rep(1, 10))
)

test_that("no covariate, no control_env returns 0.5 on df_test4", {
  gp <- suppressWarnings(
    calculate_grand_plasticity(
      trait_values = df_test4$Column0,
      env_data      = df_test4$Column1
    )
  )
  expect_equal(gp, 0.5)
})

test_that("with matching covariate yields zero plasticity", {
  gp <- suppressWarnings(
    calculate_grand_plasticity(
      trait_values  = df_test4$Column0,
      env_data       = df_test4$Column1,
      covariate_data = df_test4$Column2
    )
  )
  expect_equal(gp, 0.5)
})


test_that("excluding control_env works correctly", {
  # for df_test4, env levels are 1,2,3; control_env=2 excludes emmean=3
  # treatment emmeans = c(2,1) so sd=0.7071068, mean=1.5 => gp≈0.4714045
  gp <- suppressWarnings(calculate_grand_plasticity(
    trait_values   = df_test4$Column0,
    env_data        = df_test4$Column1,
    control_env     = 2
  ))
  expect_equal(gp, sd(c(2,1))/mean(c(2,1)))
})



