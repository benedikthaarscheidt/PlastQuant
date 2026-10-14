# =============================================================================
# simulate_reaction_norms.R — simulate SNPs and, from them, the reaction norms
#                             (mother norms and replicates) of every shape
# =============================================================================
# simulate_reaction_norms(cfg) runs the model in this order:
#   SNP genotypes and causal SNPs -> SNP scores -> genetic values of the mother
#   norms -> repeatability -> latent values of the plants -> parameter values ->
#   reaction norms -> measurements.
# The result has the same layout as own data (see pipeline_input.R), plus the
# simulation truth and the intermediate values needed by variance_report().
# =============================================================================

simulate_reaction_norms <- function(cfg) {
  # Genetics: drawn first and from their own seed, so they are identical for every h2
  set.seed(cfg$seed)
  geno <- simulate_snp_genotypes(length(cfg$shapes) * cfg$n_per_shape, cfg$n_snps)
  causal <- draw_causal_snps(cfg)

  # All standard-normal draws come from a second stream, also independent of h2
  set.seed(cfg$seed + 1)
  t <- (seq_len(cfg$n_points) - 1) / (cfg$n_points - 1)

  per_shape <- lapply(seq_along(cfg$shapes), function(k) {
    rows <- (k - 1) * cfg$n_per_shape + seq_len(cfg$n_per_shape)
    simulate_shape(cfg, cfg$shapes[k], geno, causal, rows, t)
  })
  names(per_shape) <- cfg$shapes

  list(
    snp_matrix        = geno$X,
    environments      = 1 + 9 * t,
    phenotypes        = lapply(per_shape, `[[`, "phenotypes"),
    truth             = causal[, c("Shape", "Parameter", "SNP", "Effect")],
    mother_curves     = lapply(per_shape, `[[`, "mother_curves"),
    mother_parameters = do.call(rbind, lapply(per_shape, `[[`, "mother_parameters")),
    plant_parameters  = lapply(per_shape, `[[`, "plant_parameters"),
    latent            = lapply(per_shape, `[[`, "latent"),
    measurement_sd    = vapply(per_shape, `[[`, numeric(1), "measurement_sd"),
    config            = cfg
  )
}

# One shape: genotypes in `rows` of the SNP matrix
simulate_shape <- function(cfg, shape, geno, causal, rows, t) {
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

  # Reaction norms at the measurement points
  mother_Y <- reaction_norm_matrix(shape, t, lapply(mother, as.matrix), cfg)
  plant_Y <- reaction_norm_matrix(shape, t, plant, cfg)

  # Measurement
  sigma_e <- measurement_error_sd(plant_Y, cfg$error_share)
  Y_obs <- add_measurement_error(plant_Y, sigma_e, matrix(rnorm(length(plant_Y)), nrow(plant_Y)))

  plant_id <- rep(ids, times = cfg$n_rep)
  plant_rep <- rep(seq_len(cfg$n_rep), each = n)
  ord <- order(match(plant_id, ids), plant_rep)
  env <- 1 + 9 * t

  list(
    phenotypes = long_table(plant_id[ord], plant_rep[ord], env, Y_obs[ord, , drop = FALSE]),
    mother_curves = long_table(ids, rep(0L, n), env, mother_Y)[, c("Genotype", "Environment", "Trait")],
    mother_parameters = parameter_table(data.frame(Genotype = ids, Shape = shape), mother, rho, cfg),
    plant_parameters = data.frame(
      Genotype = plant_id[ord], Replicate = plant_rep[ord],
      as.data.frame(lapply(plant, function(m) as.vector(m)[ord])), rho = rep(rho, times = cfg$n_rep)[ord]
    ),
    latent = list(snp_score = S, genetic_value = ell, tau = tau),
    measurement_sd = sigma_e
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
