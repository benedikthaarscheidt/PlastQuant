# =============================================================================
# genetic_values.R — genetic value of a genotype, its repeatability and the
#                    latent values of its replicate plants
# =============================================================================

# ell_i = sqrt(h2) * S_i + sqrt(1 - h2) * B_i  (variance 1 for every h2)
genetic_value <- function(snp_score, background, h2) sqrt(h2) * snp_score + sqrt(1 - h2) * background

# Percentile of a latent value; clamped so that extreme values never give exactly 0 or 1
latent_to_unit <- function(ell) {
  u <- pnorm(ell)
  u[] <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  u
}

# rho_i = Q_Beta(ell_iV): correlation between two replicates of genotype i
repeatability_from_latent <- function(ell_V, beta_shape) {
  qbeta(latent_to_unit(ell_V), beta_shape[1], beta_shape[2])
}

# tau_i = sqrt((1 - rho_i) / rho_i): spread of the replicates around the mother norm
replicate_spread <- function(rho) sqrt((1 - rho) / rho)

# ell_ir = ell_i + tau_i * zeta_ir; zeta is a genotypes x replicates matrix of N(0, 1)
plant_latent_values <- function(ell, tau, zeta) ell + tau * zeta
