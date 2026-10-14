# =============================================================================
# variance_report.R — where the variance of a simulation comes from
# =============================================================================
# Two levels:
#   latent: per parameter, the level on which h2 and the repeatability are set
#           (SNP share of the genetic variance between mother norms = h2)
#   trait:  per measurement y(t), the variance that results from the model:
#           causal SNPs, genetic background, their non-linear interaction,
#           replicates within genotypes and measurement error
# =============================================================================

variance_report <- function(sim, n_mc = 2000) {
  cfg <- sim$config
  latent <- do.call(rbind, lapply(cfg$shapes, function(s) {
    lat <- sim$latent[[s]]
    do.call(rbind, lapply(reaction_norm_parameters(cfg, s), function(theta) {
      data.frame(
        Shape = s, Parameter = theta, target_h2 = cfg$h2,
        snp_share = var(sqrt(cfg$h2) * lat$snp_score[[theta]]) / var(lat$genetic_value[[theta]]),
        mean_repeatability = mean(1 / (1 + lat$tau^2)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  trait <- do.call(rbind, lapply(cfg$shapes, function(s) {
    d <- trait_variance_decomposition(cfg, s, n_mc)
    me <- sim$measurement_sd[[s]]^2
    total <- d$total + me
    data.frame(
      Shape = s,
      Component = c("causal SNPs", "genetic background", "SNP x background", "within genotype", "measurement error"),
      Share = c(d$snp, d$background, d$interaction, d$within, me) / total,
      stringsAsFactors = FALSE
    )
  }))
  measurement <- data.frame(
    Shape = cfg$shapes,
    sigma_e = sim$measurement_sd[cfg$shapes],
    n_negative = vapply(cfg$shapes, function(s) sum(sim$phenotypes[[s]]$Trait < 0), integer(1)),
    row.names = NULL
  )
  list(latent = latent, trait = trait, measurement = measurement)
}

# Monte Carlo (pick-freeze) decomposition of Var(y(t)) averaged over the measurement
# points, for a large population drawn from the same model (n_mc genotypes).
trait_variance_decomposition <- function(cfg, shape, n_mc = 2000) {
  t <- (seq_len(cfg$n_points) - 1) / (cfg$n_points - 1)
  traits <- c(reaction_norm_parameters(cfg, shape), "V")
  draw <- function() setNames(lapply(traits, function(x) rnorm(n_mc)), traits)
  curves <- function(S, B, Z) {
    ell <- lapply(traits, function(x) genetic_value(S[[x]], B[[x]], cfg$h2))
    names(ell) <- traits
    tau <- replicate_spread(repeatability_from_latent(ell$V, cfg$repeatability_beta))
    p <- lapply(setdiff(traits, "V"), function(x) {
      latent_to_parameter(plant_latent_values(ell[[x]], tau, Z[[x]]), x, shape, cfg)
    })
    names(p) <- setdiff(traits, "V")
    reaction_norm_matrix(shape, t, p, cfg)
  }
  S <- draw(); B <- draw(); S2 <- draw(); B2 <- draw(); Z1 <- draw(); Z2 <- draw()
  y <- curves(S, B, Z1)
  same_genotype <- curves(S, B, Z2)
  same_snps <- curves(S, B2, Z2)
  same_background <- curves(S2, B, Z2)
  avg_cov <- function(a, b) mean(vapply(seq_along(t), function(j) cov(a[, j], b[, j]), numeric(1)))
  total <- mean(apply(y, 2, var))
  genetic <- avg_cov(y, same_genotype)
  snp <- avg_cov(y, same_snps)
  background <- avg_cov(y, same_background)
  list(total = total, snp = snp, background = background,
       interaction = genetic - snp - background, within = total - genetic)
}
