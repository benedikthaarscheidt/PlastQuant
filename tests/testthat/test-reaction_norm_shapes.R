trapz01 <- function(f) (sum(f) - (f[1] + f[length(f)]) / 2) / (length(f) - 1)
shapes <- list(linear = list(), gaussian = list(c = 0.3, w = 0.08),
               sinusoidal = list(k = 1.4, psi = 0.7), wave = list(k = 1.6, c = 0.5, w = 0.4))

test_that("normalised shapes have mean 0 and largest deviation 1", {
  for (s in names(shapes)) {
    v <- normalised_shape(s, normalisation_grid(), shapes[[s]])
    expect_equal(trapz01(v), 0, tolerance = 1e-10, info = s)
    expect_equal(max(abs(v)), 1, tolerance = 1e-12, info = s)
  }
})

test_that("linear normalised shape is 2t - 1", {
  t <- c(0, 0.2, 0.5, 1)
  expect_equal(normalised_shape("linear", t, list()), 2 * t - 1)
})

test_that("L is the curve mean, a the largest relative deviation, y > 0", {
  y <- reaction_norm_curve("gaussian", normalisation_grid(), L = 1.3, a = 0.4, g = shapes$gaussian)
  expect_equal(trapz01(y), 1.3, tolerance = 1e-10)
  expect_equal(max(abs(y - 1.3)) / 1.3, 0.4, tolerance = 1e-12)
  expect_true(min(y) >= 1.3 * 0.6 - 1e-12)
})

test_that("unknown shape is an error", {
  expect_error(shape_curve("cubic", 0.5, list()), "unknown shape")
})
