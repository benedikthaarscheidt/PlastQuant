# =============================================================================
# reaction_norm_shapes.R — raw shapes, their normalisation and the reaction norm
# =============================================================================

# Fine grid on [0, 1] for the normalisation; its step 1/2000 contains the
# measurement points t_j = (j-1)/(T-1) for T in {2, 3, 5, 6, 9, 11, ...}
normalisation_grid <- function() seq(0, 1, length.out = 2001)

shape_curve <- function(shape, t, g) {
  switch(shape,
    linear     = t,
    gaussian   = exp(-(t - g$c)^2 / (2 * g$w^2)),
    sinusoidal = sin(2 * pi * g$k * t + g$psi),
    wave       = sin(2 * pi * g$k * (t - g$c)) * exp(-(t - g$c)^2 / (2 * g$w^2)),
    stop("unknown shape: ", shape)
  )
}

# Mean of the raw curve over [0, 1] (trapezoid rule) and its largest deviation from it
shape_normalisation <- function(shape, g) {
  f <- shape_curve(shape, normalisation_grid(), g)
  n <- length(f)
  mean_f <- (sum(f) - (f[1] + f[n]) / 2) / (n - 1)
  dmax <- max(abs(f - mean_f))
  if (!(dmax > 0)) stop("shape is constant on [0, 1]: ", shape)
  list(mean = mean_f, dmax = dmax)
}

# F°(t): mean 0 and largest deviation 1, so it carries only the shape
normalised_shape <- function(shape, t, g) {
  k <- shape_normalisation(shape, g)
  (shape_curve(shape, t, g) - k$mean) / k$dmax
}

# y(t) = L (1 + a F°(t)): L is the mean of the curve, a its largest relative deviation
reaction_norm_curve <- function(shape, t, L, a, g) L * (1 + a * normalised_shape(shape, t, g))
