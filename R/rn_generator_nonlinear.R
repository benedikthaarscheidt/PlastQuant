#' Simulate Nonlinear Reaction Norm Data for Multiple Patterns
#'
#' Generates synthetic nonlinear reaction norm data (Gaussian, Sinusoidal, Wave) for multiple
#' genotypes and replicates, optionally injects outliers, plots examples, and saves both
#' data and figures to disk.
#'
#' @param num_genotypes Integer: number of genotypes to simulate (default 20).
#' @param replicates Integer: number of replicates per genotype (default 3).
#' @param env_min Numeric: minimum environment value (default 1).
#' @param env_max Numeric: maximum environment value (default 10).
#' @param env_length Integer: number of points along the environmental gradient (default 50).
#' @param var_offset Named numeric vector of length 3: variance of intercept random effect
#'        for each norm type (default c(gaussian=0.5, sinusoidal=0.3, wave=0.7)).
#' @param var_slope Named numeric vector of length 3: variance of slope random effect
#'        for each norm type (default c(gaussian=0.1, sinusoidal=0.05, wave=0.15)).
#' @param correlation Named numeric vector of length 3: correlation of random intercept and slope
#'        for each norm type (default c(gaussian=-0.5, sinusoidal=-0.3, wave=-0.7)).
#' @param outliers Character: one of "none","positive","negative","both"
#'        to inject a single outlier into the mother curve (default "none").
#' @param seed Integer: random seed for reproducibility (default 12).
#' @param plot Logical: whether to generate example plots (default TRUE).
#' @param genotypes_per_plot Integer: number of genotypes per example plot (default 3).
#' @param plot_file Character or NULL: filepath for combined PDF of plots. If NULL, no file written.
#' @param csv_dir Character or NULL: directory path to save CSV files. If NULL, no files written.
#'
#' @return A data.frame with all simulated observations (Genotype, Replicate,
#'         Environment, Trait, ReactionNorm, BaseShift, Slope, VarianceOffset,
#'         VarianceSlope, Covariance) is invisibly returned.
#'
#' @details
#' - Three pattern functions are used:
#'   * Gaussian: amplitude * exp(-width*(env-center)^2)
#'   * Sinusoidal: amplitude * sin(frequency*env + phase)
#'   * Wave: amplitude * sin(frequency*(env-5)) * exp(-0.1*(env-5)^2)
#' - Random intercept/slope deviations are drawn via MASS::mvrnorm.
#' - Mother curves may receive a single outlier modification.
#' - Example plots are created via ggplot2 and optionally saved.
#' - CSVs for each pattern and combined data are optionally saved.
#'
#' @importFrom MASS mvrnorm
#' @import ggplot2 dplyr
#' @export
simulate_nonlin_reaction_norms <- function(
    num_genotypes      = 20,
    replicates         = 3,
    env_min            = 1,
    env_max            = 10,
    env_length         = 50,
    var_offset         = c(gaussian=0.5, sinusoidal=0.3, wave=0.7),
    var_slope          = c(gaussian=0.1, sinusoidal=0.05, wave=0.15),
    correlation        = c(gaussian=-0.5, sinusoidal=-0.3, wave=-0.7),
    outliers           = "none",
    seed               = 12,
    plot               = TRUE,
    genotypes_per_plot = 3,
    plot_file          = NULL,
    csv_dir            = NULL
) {
  # Helpers -------------------------------------------------------------------
  build_cov_matrix <- function(vo, vs, cor) {
    cov <- cor * sqrt(vo * vs)
    matrix(c(vo, cov, cov, vs), nrow = 2)
  }
  generate_gaussian <- function(env, bs, sl) bs + sl * 4 * exp(-0.2*(env-5)^2)
  generate_sinus   <- function(env, bs, sl) bs + sl * 4 * sin(0.5*env)
  generate_wave    <- function(env, bs, sl) bs + sl * 7 * sin(1.2*(env-5)) * exp(-0.1*(env-5)^2)
  inject_outlier   <- function(curve) {
    i <- sample(seq_along(curve), 1)
    factor <- switch(outliers,
                     positive=runif(1,1,2),
                     negative=runif(1,0,1),
                     both    =runif(1,0,4),
                     1)
    curve[i] <- curve[i] * factor
    curve
  }
  ensure_dir <- function(path) {
    d <- dirname(path)
    if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  }

  # Setup ---------------------------------------------------------------------
  set.seed(seed)
  env <- seq(env_min, env_max, length.out=env_length)
  base_shifts <- runif(num_genotypes, 10, 30)
  slopes      <- runif(num_genotypes, 0.5, 1.5) * sample(c(-1,1), num_genotypes, TRUE)

  # Generate per-genotype data -----------------------------------------------
  gen_data <- function(norm_fun, label) {
    do.call(rbind, lapply(seq_len(num_genotypes), function(id) {
      bs <- base_shifts[id]; sl <- slopes[id]
      covm <- build_cov_matrix(var_offset[label], var_slope[label], correlation[label])
      mother <- norm_fun(env, bs, sl)
      if (outliers != "none") mother <- inject_outlier(mother)
      base   <- norm_fun(env, bs, sl)
      mom_df <- data.frame(
        Genotype     = id,
        Replicate    = 0,
        Environment  = env,
        Trait        = mother,
        ReactionNorm = label,
        BaseShift    = bs,
        Slope        = sl,
        VarianceOffset = NA,
        VarianceSlope  = NA,
        Covariance     = NA
      )
      reps <- do.call(rbind, lapply(seq_len(replicates), function(rp) {
        re <- MASS::mvrnorm(1, mu=c(0,0), Sigma=covm)
        pc <- norm_fun(env, bs+re[1], sl+re[2])
        dev <- pc - base
        rc  <- mother + dev
        data.frame(
          Genotype     = id,
          Replicate    = rp,
          Environment  = env,
          Trait        = rc,
          ReactionNorm = label,
          BaseShift    = bs + re[1],
          Slope        = sl + re[2],
          VarianceOffset = var_offset[label],
          VarianceSlope  = var_slope[label],
          Covariance     = correlation[label]*sqrt(var_offset[label]*var_slope[label])
        )
      }))
      rbind(mom_df, reps)
    }))
  }

  gaussian_df <- gen_data(generate_gaussian,  "gaussian")
  sinusoidal_df <- gen_data(generate_sinus,   "sinusoidal")
  wave_df <- gen_data(generate_wave,          "wave")
  all_df <- rbind(gaussian_df, sinusoidal_df, wave_df)

  # Plot ----------------------------------------------------------------------
  if (plot && !is.null(plot_file)) {
    # Ensure output directory exists
    ensure_dir(plot_file)
    # Open multi-page PDF device
    grDevices::pdf(plot_file, width = 8, height = 6)

    # Choose subset of genotypes
    chosen <- sample(unique(all_df$Genotype), genotypes_per_plot)

    # Plot function with dashed replicates
    plot_panel <- function(df, title) {
      df_sub <- df[df$Genotype %in% chosen, ]
      ggplot2::ggplot(df_sub, ggplot2::aes(
        x = Environment,
        y = Trait,
        color = factor(Genotype),
        group = interaction(Genotype, Replicate),
        linetype = Replicate == 0,
        alpha = Replicate == 0
      )) +
        ggplot2::geom_line(size = 1) +
        ggplot2::scale_linetype_manual(
          values = c("TRUE" = "solid", "FALSE" = "dashed"),
          guide = "none"
        ) +
        ggplot2::scale_alpha_manual(
          values = c("TRUE" = 1, "FALSE" = 0.6),
          guide = "none"
        ) +
        ggplot2::labs(title = title, x = "Environment", y = "Trait") +
        ggplot2::theme_minimal()
    }

    # Gaussian Norms
    print(plot_panel(gaussian_df,  "Gaussian Norms"))
    # Sinusoidal Norms
    print(plot_panel(sinusoidal_df, "Sinusoidal Norms"))
    # Wave Norms
    print(plot_panel(wave_df,       "Wave Norms"))

    # Close device
    grDevices::dev.off()
  }

  # Write CSVs ---------------------------------------------------------------- ---------------------------------------------------------------- ----------------------------------------------------------------
  if (!is.null(csv_dir)) {
    if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)
    utils::write.csv(gaussian_df, file.path(csv_dir, "gaussian_data.csv"), row.names = FALSE)
    utils::write.csv(sinusoidal_df, file.path(csv_dir, "sinusoidal_data.csv"), row.names = FALSE)
    utils::write.csv(wave_df, file.path(csv_dir, "wave_data.csv"), row.names = FALSE)
    utils::write.csv(all_df, file.path(csv_dir, "all_combined_data.csv"), row.names = FALSE)
  }

  invisible(all_df)
}
