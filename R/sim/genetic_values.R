# =============================================================================
# genetic_values.R — genetic value of a genotype, the stability of its replicates
#                    and the latent values of its replicate plants
# =============================================================================
# Latent scale of every reaction norm parameter theta (genotype i, replicate plant r):
#   G_i       = sqrt(h2) * (sqrt(q) * S_i + sqrt(1 - q) * B_i)   genetic value, Var = V_A = h2
#   l_iV      = sqrt(q) * S_iV + sqrt(1 - q) * B_iV               genetic value of the stability V, Var = 1
#   sigma_i^2 = exp(s * l_iV - s^2 / 2)                            relative spread of the replicates, E = 1
#   l_ir      = G_i + sqrt((1 - h2) * sigma_i^2) * zeta_ir         plant, Var = V_P = 1
# S: SNP score of the causal SNPs, B: polygenic background (additive), zeta: deviation
# of the plant; all N(0, 1) and independent. V_W = (1 - h2) E[sigma_i^2] = 1 - h2, so
# h2 = V_A / V_P is the narrow-sense heritability for every s, and q * h2 the share of
# the causal SNPs (q * h2 / Mc per causal SNP). The stability SNPs (V) change only the
# spread of the replicates, not their mean, and add nothing to V_A.
# =============================================================================

# G_i = sqrt(h2) * (sqrt(q) * S_i + sqrt(1 - q) * B_i); h2 = 1 gives the stability value l_iV
genetic_value <- function(snp_score, background, h2, q) sqrt(h2) * (sqrt(q) * snp_score + sqrt(1 - q) * background)

# Percentile of a latent value; clamped so that extreme values never give exactly 0 or 1
latent_to_unit <- function(ell) {
  u <- pnorm(ell)
  u[] <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  u
}

# sigma_i^2 = exp(s * l_iV - s^2 / 2): spread of the replicates of genotype i relative to
# the average (mean 1 over the genotypes); s = 0: all genotypes equally stable
replicate_spread <- function(ell_V, s) exp(s * ell_V - s^2 / 2)

# rho_i = h2 / (h2 + (1 - h2) sigma_i^2): correlation of two replicates at stability sigma_i^2
repeatability <- function(sigma2, h2) h2 / (h2 + (1 - h2) * sigma2)

# l_ir = G_i + sqrt((1 - h2) sigma_i^2) * zeta_ir; zeta is a genotypes x replicates matrix
# of N(0, 1) (or a vector, one plant per genotype)
plant_latent_values <- function(G, sigma2, zeta, h2) G + sqrt((1 - h2) * sigma2) * zeta

# Settings saved before this model (no cfg$q: h2 was the share of the causal SNPs in the
# genetic variance 1, and rho_i ~ Beta(repeatability_beta)); kept so that populations
# simulated with them can still be measured at new resolutions (variance_report.R)
legacy_plant_latent_values <- function(S, B, Z, cfg) {
  ell <- lapply(setNames(names(S), names(S)), function(x) sqrt(cfg$h2) * S[[x]] + sqrt(1 - cfg$h2) * B[[x]])
  rho <- qbeta(latent_to_unit(ell$V), cfg$repeatability_beta[1], cfg$repeatability_beta[2])
  tau <- sqrt((1 - rho) / rho)
  thetas <- setdiff(names(S), "V")
  lapply(setNames(thetas, thetas), function(x) ell[[x]] + tau * Z[[x]])
}
