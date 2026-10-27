# =============================================================================
# variance_report.R — where the variance of a simulation comes from
# =============================================================================
# Two levels:
#   latent: per parameter, the level on which the heritability is set: target and
#           realised V_A, V_Q (causal SNPs), V_W (replicates), h2 = V_A / V_P and the
#           share of the causal SNPs (see genetic_values.R)
#   trait:  per measurement y(t), the variance that results from the model:
#           causal SNPs, genetic background, their non-linear interaction,
#           replicates within genotypes and measurement error
# =============================================================================

variance_report <- function(sim, n_mc = 2000) {
  c(list(latent = latent_variance(sim)), trait_variance(sim, n_mc))
}

# Latent level: needs only the reaction norms (simulate_norms() or a measured simulation).
# Targets (exact in the population) next to the values realised over the genotypes of
# the population: V_A = Var(G), V_Q = Var(sqrt(q h2) S), V_W = (1 - h2) mean(sigma_i^2).
latent_variance <- function(norms) {
  cfg <- norms$config
  do.call(rbind, lapply(cfg$shapes, function(s) {
    lat <- norms$latent[[s]]
    V_W <- (1 - cfg$h2) * mean(lat$sigma2)
    do.call(rbind, lapply(reaction_norm_parameters(cfg, s), function(theta) {
      V_A <- var(lat$genetic_value[[theta]])
      V_Q <- var(sqrt(cfg$q * cfg$h2) * lat$snp_score[[theta]])
      data.frame(
        Shape = s, Parameter = theta,
        target_h2 = cfg$h2, h2 = V_A / (V_A + V_W),
        target_snp_share = cfg$q * cfg$h2, snp_share = V_Q / (V_A + V_W),
        target_per_snp = cfg$q * cfg$h2 / cfg$Mc,
        V_A = V_A, V_Q = V_Q, V_W = V_W,
        mean_repeatability = mean(repeatability(lat$sigma2, cfg$h2)),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# Trait level of one measurement of the norms (measure_norms()): the decomposition at its
# sampling points, the measurement error and the number of negative values
trait_variance <- function(measured, n_mc = 2000) {
  cfg <- measured$config
  t <- (measured$environments - 1) / 9
  trait <- do.call(rbind, lapply(cfg$shapes, function(s) {
    d <- trait_variance_decomposition(cfg, s, n_mc, t)
    me <- measured$measurement_sd[[s]]^2
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
    sigma_e = unname(measured$measurement_sd[cfg$shapes]),
    n_negative = vapply(cfg$shapes, function(s) sum(measured$phenotypes[[s]]$Trait < 0), integer(1)),
    row.names = NULL
  )
  list(trait = trait, measurement = measurement)
}

# Monte Carlo (pick-freeze) decomposition of Var(y(t)) averaged over the measurement
# points t (default: cfg$n_points equidistant points), for a large population drawn
# from the same model (n_mc genotypes).
trait_variance_decomposition <- function(cfg, shape, n_mc = 2000, t = (seq_len(cfg$n_points) - 1) / (cfg$n_points - 1)) {
  traits <- c(reaction_norm_parameters(cfg, shape), "V")
  draw <- function() setNames(lapply(traits, function(x) rnorm(n_mc)), traits)
  thetas <- setdiff(traits, "V")
  curves <- function(S, B, Z) {
    plant <- if (is.null(cfg$q)) {
      legacy_plant_latent_values(S, B, Z, cfg)   # settings saved before the narrow-sense model
    } else {
      sigma2 <- replicate_spread(genetic_value(S$V, B$V, 1, cfg$q), cfg$stability_sd)
      lapply(setNames(thetas, thetas), function(x) {
        plant_latent_values(genetic_value(S[[x]], B[[x]], cfg$h2, cfg$q), sigma2, Z[[x]], cfg$h2)
      })
    }
    p <- lapply(setNames(thetas, thetas), function(x) latent_to_parameter(plant[[x]], x, shape, cfg))
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
