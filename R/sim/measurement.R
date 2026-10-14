# =============================================================================
# measurement.R — additive measurement error
# =============================================================================

# SD of the error so that it makes up `share` of the trait variance,
# averaged over the measurement points (Y: plants x measurement points)
measurement_error_sd <- function(Y, share) sqrt(share / (1 - share) * mean(apply(Y, 2, var)))

# Y_obs = Y + sigma_e * noise; noise is a matrix of N(0, 1) of the same size
add_measurement_error <- function(Y, sigma_e, noise) Y + sigma_e * noise
