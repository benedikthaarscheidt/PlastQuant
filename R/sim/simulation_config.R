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

# Shape parameters of each reaction-norm shape (besides level L and plasticity a):
#   c    centre of the peak on [0, 1]
#   w    width of the peak
#   k    frequency (cycles over [0, 1])
#   psi  phase of the sine
shape_parameters <- function() {
  list(
    linear     = character(0),        # only L and a
    gaussian   = c("c", "w"),         # peak at c with width w
    sinusoidal = c("k", "psi"),       # sine with frequency k and phase psi
    wave       = c("k", "c", "w")     # sine of frequency k under a Gaussian envelope at c, width w
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

sim_config <- function(h2,                                  # share of the latent genetic variance from the causal SNPs (0..1)
                       Mc,                                  # causal SNPs per parameter (and for V) of every shape
                       n_per_shape = 200,                   # genotypes (lines) per shape
                       n_rep = 3,                           # replicate plants per genotype
                       n_snps = 200000,                     # simulated SNPs (causal + neutral)
                       n_points = 6,                        # measurement points per reaction norm; environments 1 + 9 t
                       seed = 1,                            # seed for the genetics; seed + 1 for all other normal draws
                       sigma_L = 0.1,                       # SD of ln(L): spread of the level, L = exp(sigma_L * ell)
                       repeatability_beta = c(18, 2),       # Beta shape of the repeatability rho over genotypes (mean 0.9)
                       error_share = 0.05,                  # measurement error as share of the trait variance
                       study_space = default_study_space()  # range and distribution of every bounded parameter
                       ) {
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
reaction_norm_parameters <- function(cfg, shape) {
  if (is.null(cfg$shape_parameters)) stop("cfg has no shape_parameters; build the settings with sim_config()")
  if (!shape %in% names(cfg$shape_parameters)) stop("unknown shape: ", shape)
  c("L", "a", cfg$shape_parameters[[shape]])
}
