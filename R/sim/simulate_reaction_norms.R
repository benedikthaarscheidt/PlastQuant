# =============================================================================
# simulate_reaction_norms.R — simulate SNPs and, from them, the reaction norms
#                             (mother norms and replicates) of every shape
# =============================================================================
# Two separate steps:
#   simulate_norms(cfg)       generates the population: SNP genotypes and causal SNPs ->
#                             SNP scores -> genetic values of the mother norms ->
#                             repeatability -> latent values of the plants -> parameter
#                             values of every mother norm and replicate plant. A reaction
#                             norm is fully given by its parameters.
#   measure_norms(norms, n)   samples these same reaction norms at n equidistant points
#                             and adds the measurement error.
# simulate_reaction_norms(cfg) does both, at cfg$n_points (the population's resolution).
# Its result has the same layout as own data (see pipeline_input.R), plus the
# simulation truth and the intermediate values needed by variance_report().
# =============================================================================

simulate_reaction_norms <- function(cfg) {
  measure_norms(simulate_norms(cfg), cfg$n_points)
}

simulate_norms <- function(cfg) {
  # Genetics: drawn first and from their own seed, so they are identical for every h2
  set.seed(cfg$seed)
  geno <- simulate_snp_genotypes(length(cfg$shapes) * cfg$n_per_shape, cfg$n_snps)
  causal <- draw_causal_snps(cfg)

  # All standard-normal draws of the norms come from a second stream, also independent of h2
  set.seed(cfg$seed + 1)
  per_shape <- lapply(seq_along(cfg$shapes), function(k) {
    rows <- (k - 1) * cfg$n_per_shape + seq_len(cfg$n_per_shape)
    simulate_shape_norms(cfg, cfg$shapes[k], geno, causal, rows)
  })
  names(per_shape) <- cfg$shapes

  list(
    snp_matrix        = geno$X,
    truth             = causal[, c("Shape", "Parameter", "SNP", "Effect")],
    mother_parameters = do.call(rbind, lapply(per_shape, `[[`, "mother_parameters")),
    plant_parameters  = lapply(per_shape, `[[`, "plant_parameters"),
    latent            = lapply(per_shape, `[[`, "latent"),
    config            = cfg
  )
}

# Sample the reaction norms of simulate_norms() at n_points equidistant points
# t_j = (j - 1) / (n_points - 1), environments 1 + 9 t, and add measurement error with
# error share cfg$error_share of the trait variance at these points (error = FALSE:
# the noise-free values). The error comes from its own stream (seed + 2), so the same
# design always gives the same measurements.
measure_norms <- function(norms, n_points, error = TRUE) {
  cfg <- norms$config
  if (n_points < 2) stop("n_points must be at least 2")
  t <- (seq_len(n_points) - 1) / (n_points - 1)
  env <- 1 + 9 * t
  set.seed(cfg$seed + 2)
  per_shape <- lapply(cfg$shapes, function(shape) {
    thetas <- reaction_norm_parameters(cfg, shape)
    pp <- norms$plant_parameters[[shape]]
    mp <- norms$mother_parameters[norms$mother_parameters$Shape == shape, ]
    plant_Y <- reaction_norm_matrix(shape, t, as.list(pp[thetas]), cfg)
    mother_Y <- reaction_norm_matrix(shape, t, as.list(mp[thetas]), cfg)
    sigma_e <- measurement_error_sd(plant_Y, cfg$error_share)
    noise <- matrix(rnorm(length(plant_Y)), nrow(plant_Y))
    Y_obs <- if (error) add_measurement_error(plant_Y, sigma_e, noise) else plant_Y
    list(
      phenotypes = long_table(pp$Genotype, pp$Replicate, env, Y_obs),
      mother_curves = long_table(mp$Genotype, rep(0L, nrow(mp)), env, mother_Y)[, c("Genotype", "Environment", "Trait")],
      measurement_sd = sigma_e
    )
  })
  names(per_shape) <- cfg$shapes
  c(norms, list(
    environments   = env,
    phenotypes     = lapply(per_shape, `[[`, "phenotypes"),
    mother_curves  = lapply(per_shape, `[[`, "mother_curves"),
    measurement_sd = vapply(per_shape, `[[`, numeric(1), "measurement_sd")
  ))
}

# Parameters of the mother norms and replicate plants of one shape: genotypes in
# `rows` of the SNP matrix
simulate_shape_norms <- function(cfg, shape, geno, causal, rows) {
  thetas <- reaction_norm_parameters(cfg, shape)
  ids <- rownames(geno$X)[rows]
  n <- length(rows)

  # Genetic value of the mother norm for every parameter and for V
  S <- list(); ell <- list()
  for (theta in c(thetas, "V")) {
    cs <- causal[causal$Shape == shape & causal$Parameter == theta, ]
    S[[theta]] <- snp_score(geno$X, geno$p, rows, cs$SNP_index, cs$Effect)
    ell[[theta]] <- genetic_value(S[[theta]], rnorm(n), cfg$h2)
  }

  # Repeatability and spread of the replicates
  rho <- repeatability_from_latent(ell$V, cfg$repeatability_beta)
  tau <- replicate_spread(rho)

  # Latent values of the plants and parameter values (genotypes x replicates)
  plant <- list(); mother <- list()
  for (theta in thetas) {
    zeta <- matrix(rnorm(n * cfg$n_rep), n, cfg$n_rep)
    plant[[theta]] <- latent_to_parameter(plant_latent_values(ell[[theta]], tau, zeta), theta, shape, cfg)
    mother[[theta]] <- latent_to_parameter(ell[[theta]], theta, shape, cfg)
  }

  # one row per plant, genotype-major, replicates in order
  plant_id <- rep(ids, times = cfg$n_rep)
  plant_rep <- rep(seq_len(cfg$n_rep), each = n)
  ord <- order(match(plant_id, ids), plant_rep)

  list(
    mother_parameters = parameter_table(data.frame(Genotype = ids, Shape = shape), mother, rho, cfg),
    plant_parameters = data.frame(
      Genotype = plant_id[ord], Replicate = plant_rep[ord],
      as.data.frame(lapply(plant, function(m) as.vector(m)[ord])), rho = rep(rho, times = cfg$n_rep)[ord],
      stringsAsFactors = FALSE
    ),
    latent = list(snp_score = S, genetic_value = ell, tau = tau)
  )
}

# Curves of all plants (rows: genotype-major within each replicate column of the
# parameter matrices) at the points t
reaction_norm_matrix <- function(shape, t, params, cfg) {
  n_plants <- length(params$L)
  shape_par <- cfg$shape_parameters[[shape]]
  Y <- matrix(NA_real_, n_plants, length(t))
  for (j in seq_len(n_plants)) {
    g <- lapply(shape_par, function(theta) params[[theta]][j])
    names(g) <- shape_par
    Y[j, ] <- reaction_norm_curve(shape, t, params$L[j], params$a[j], g)
  }
  Y
}

# Long table: one row per plant and measurement point
long_table <- function(genotype, replicate, env, Y) {
  data.frame(
    Genotype    = rep(genotype, each = length(env)),
    Replicate   = rep(as.integer(replicate), each = length(env)),
    Environment = rep(env, times = length(genotype)),
    Trait       = as.vector(t(Y)),
    stringsAsFactors = FALSE
  )
}

# Genotype table with one column per parameter of any shape (NA where the shape
# has no such parameter), so that the tables of all shapes can be stacked
parameter_table <- function(df, values, rho, cfg) {
  all_params <- unique(c("L", "a", unlist(cfg$shape_parameters)))
  for (theta in all_params) df[[theta]] <- if (theta %in% names(values)) as.vector(values[[theta]]) else NA_real_
  df$rho <- rho
  df
}
