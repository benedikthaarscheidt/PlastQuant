# =============================================================================
# simulation_config.R — every setting of the simulation in one place
# =============================================================================
# All variance components are controlled here:
#   h2                  share of the genetic variance between mother norms that
#                       comes from the causal SNPs (the rest is genetic background)
#   repeatability_beta  Beta shape of the repeatability rho over genotypes; it sets
#                       how strongly replicates scatter around their mother norm
#   sigma_L             SD of ln(L) over mother norms (level)
#   error_share         measurement error as share of the trait variance
#   study_space         range and distribution of every bounded parameter
# =============================================================================

# Shape parameters of each reaction-norm shape (besides level L and plasticity a)
shape_parameters <- function() {
  list(
    linear     = character(0),
    gaussian   = c("c", "w"),
    sinusoidal = c("k", "psi"),
    wave       = c("k", "c", "w")
  )
}

# Range and distribution of every bounded parameter, per shape
default_study_space <- function() {
  a   <- list(type = "uniform", min = 0, max = 0.7)
  c   <- list(type = "uniform", min = 0, max = 1)
  w   <- list(type = "loguniform", min = 0.05, max = 1)
  psi <- list(type = "uniform", min = 0, max = pi)
  list(
    linear     = list(a = a),
    gaussian   = list(a = a, c = c, w = w),
    sinusoidal = list(a = a, k = list(type = "loguniform", min = 0.12, max = 2), psi = psi),
    wave       = list(a = a, k = list(type = "loguniform", min = 1, max = 2), c = c, w = w)
  )
}

sim_config <- function(h2, Mc,
                       n_per_shape = 200,
                       n_rep = 3,
                       n_snps = 200000,
                       n_points = 6,
                       seed = 1,
                       sigma_L = 0.1,
                       repeatability_beta = c(18, 2),
                       error_share = 0.05,
                       study_space = default_study_space()) {
  if (!(is.numeric(h2) && length(h2) == 1 && h2 >= 0 && h2 <= 1)) stop("h2 must be in [0, 1]")
  if (Mc < 1) stop("Mc must be at least 1")
  if (n_rep < 1) stop("n_rep must be at least 1")
  if (n_points < 2) stop("n_points must be at least 2")
  if (!(error_share >= 0 && error_share < 1)) stop("error_share must be in [0, 1)")
  shapes <- names(shape_parameters())
  list(
    h2 = h2, Mc = Mc, n_per_shape = n_per_shape, n_rep = n_rep, n_snps = n_snps,
    n_points = n_points, seed = seed, sigma_L = sigma_L,
    repeatability_beta = repeatability_beta, error_share = error_share,
    study_space = study_space, shapes = shapes, shape_parameters = shape_parameters()
  )
}

# All parameters of a reaction norm of the given shape: level, plasticity, shape
reaction_norm_parameters <- function(cfg, shape) c("L", "a", cfg$shape_parameters[[shape]])
