# =============================================================================
# simulation_config.R — every setting of the simulation in one place
# =============================================================================
# All variance components are controlled here:
#   h2                  narrow-sense heritability of every parameter on the latent
#                       scale: h2 = V_A / V_P with V_P = 1 (genetic_values.R)
#   q                   share of the additive variance V_A from the causal SNPs; the
#                       rest is polygenic background. Each causal SNP explains q h2 / Mc
#   stability_sd        how much the genotypes differ in the spread of their replicates
#                       (SD of ln sigma_i^2); changes the repeatability of single
#                       genotypes, not h2
#   sigma_L             SD of ln(L) over the plants (level)
#   error_share         measurement error as share of the trait variance
#   study_space         range and distribution of every bounded parameter
#   ld_strength         linkage disequilibrium between neighbouring SNPs (0 = none)
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

# Every argument with its default (h2 and Mc have none), to copy and adapt:
#   cfg <- sim_config(
#     h2           = 0.5,                    # narrow-sense heritability (required)
#     Mc           = 5,                      # causal SNPs per parameter (required)
#     q            = 0.5,
#     stability_sd = 0.2,
#     n_per_shape  = 200,
#     n_rep        = 3,
#     n_snps       = 200000,
#     ld_strength  = 0,
#     ld_r2        = 0.1,
#     n_points     = 6,
#     seed         = 1,
#     sigma_L      = 0.1,
#     error_share  = 0.02,
#     study_space  = default_study_space()
#   )
sim_config <- function(h2,                                  # narrow-sense heritability on the latent scale (0..1)
                       Mc,                                  # causal SNPs per parameter (and for V) of every shape
                       q = 0.5,                             # share of the additive variance from the causal SNPs (0..1)
                       stability_sd = 0.2,                  # SD of ln sigma_i^2: differences in replicate stability between genotypes (0 = none)
                       n_per_shape = 200,                   # genotypes (lines) per shape
                       n_rep = 3,                           # replicate plants per genotype
                       n_snps = 200000,                     # simulated SNPs (causal + neutral)
                       ld_strength = 0,                     # LD: probability that a SNP copies the allele of the preceding SNP on the same haplotype (0 = independent SNPs)
                       ld_r2 = 0.1,                         # r^2 that bounds the LD window (SNPs within it of a causal SNP count as finding it)
                       n_points = 6,                        # measurement points per reaction norm; environments 1 + 9 t
                       seed = 1,                            # seed for the genetics; seed + 1 for all other normal draws
                       sigma_L = 0.1,                       # SD of ln(L) over the plants, L = exp(sigma_L * ell)
                       error_share = 0.02,                  # measurement error as share of the trait variance
                       study_space = default_study_space()  # range and distribution of every bounded parameter
                       ) {
  if (!(is.numeric(h2) && length(h2) == 1 && h2 >= 0 && h2 <= 1)) stop("h2 must be in [0, 1]")
  if (!(is.numeric(q) && length(q) == 1 && q >= 0 && q <= 1)) stop("q must be in [0, 1]")
  if (!(stability_sd >= 0)) stop("stability_sd must be >= 0")
  if (Mc < 1) stop("Mc must be at least 1")
  if (n_rep < 1) stop("n_rep must be at least 1")
  if (n_points < 2) stop("n_points must be at least 2")
  if (!(error_share >= 0 && error_share < 1)) stop("error_share must be in [0, 1)")
  if (!(ld_strength >= 0 && ld_strength < 1)) stop("ld_strength must be in [0, 1)")
  if (!(ld_r2 > 0 && ld_r2 < 1)) stop("ld_r2 must be in (0, 1)")
  shapes <- names(shape_parameters())
  list(
    h2 = h2, q = q, stability_sd = stability_sd, Mc = Mc, n_per_shape = n_per_shape, n_rep = n_rep, n_snps = n_snps,
    ld_strength = ld_strength, ld_r2 = ld_r2,
    n_points = n_points, seed = seed, sigma_L = sigma_L,
    error_share = error_share,
    study_space = study_space, shapes = shapes, shape_parameters = shape_parameters()
  )
}

# All parameters of a reaction norm of the given shape: level, plasticity, shape
reaction_norm_parameters <- function(cfg, shape) {
  if (is.null(cfg$shape_parameters)) stop("cfg has no shape_parameters; build the settings with sim_config()")
  if (!shape %in% names(cfg$shape_parameters)) stop("unknown shape: ", shape)
  c("L", "a", cfg$shape_parameters[[shape]])
}
