#' Simulate and Save Synthetic Reaction Norm Data
#'
#' Generates synthetic reaction norm data for multiple genotypes and replicates,
#' with optional outlier injection, plotting, and saving to disk.
#'
#' @param num_genotypes Integer: number of genotypes to simulate (default 20).
#' @param individuals_per_genotype Integer: number of replicates per genotype (default 3).
#' @param env_min Numeric: minimum environment value (default 1).
#' @param env_max Numeric: maximum environment value (default 10).
#' @param env_length Integer: number of points along the environmental gradient (default 50).
#' @param VE Numeric: variance of intercept random effect (default 1).
#' @param VS Numeric: variance of slope random effect (default 0.2).
#' @param rES Numeric: correlation between intercept and slope random effects (default 1).
#' @param outliers Character: one of "none", "positive", "negative", or "both" for outlier injection in fixed curve (default "none").
#' @param seed Integer: random seed for reproducibility (default 42).
#' @param plot Logical: whether to generate a ggplot of a random subset of genotypes (default TRUE).
#' @param genotypes_per_plot Integer: how many genotypes to include in the example plot (default 5).
#' @param plot_file Character or NULL: filepath to save the PDF plot (if `plot = TRUE`).  If NULL, no file is written (default NULL).
#' @param csv_file Character or NULL: filepath to save the synthetic data CSV.  If NULL, no file is written (default NULL).
#'
#' @return A data.frame with columns:
#'   * Genotype, Replicate, Environment, Trait, ReactionNorm,
#'     BaseShift, Slope, VarianceOffset, VarianceSlope, Covariance
#' Invisibly returns the data.frame.
#'
#' @details
#'   - Builds a 2×2 covariance matrix from VE, VS, rES.
#'   - Generates fixed (mother) curves and adds random slope/intercept deviations for replicates.
#'   - Optionally injects a single outlier into the fixed curve.
#'   - Optionally plots a random subset of genotypes and saves as PDF.
#'   - Optionally writes the full dataset to CSV.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats lm coef
#' @import ggplot2
#' @export
simulate_linear_reaction_norms <- function(
    num_genotypes              = 20,
    individuals_per_genotype   = 3,
    env_min                    = 1,
    env_max                    = 10,
    env_length                 = 50,
    VE                         = 1,
    VS                         = 0.2,
    rES                        = 1,
    outliers                   = "none",
    seed                       = 42,
    plot                       = TRUE,
    genotypes_per_plot         = 5,
    plot_file                  = NULL,
    csv_file                   = NULL
) {
  # Helpers -------------------------------------------------------------------
  build_cov_matrix <- function(VE, VS, rES) {
    matrix(c(
      VE,        rES * sqrt(VE * VS),
      rES * sqrt(VE * VS), VS
    ), nrow = 2)
  }
  generate_fixed_norm <- function(intercept, slope, env, outliers) {
    curve <- intercept + slope * env
    curve[curve < 0] <- 0
    if (outliers == "none") return(curve)
    i <- sample(seq_along(env), 1)
    scale_factor <- switch(
      outliers,
      positive = runif(1, 1, 2),
      negative = runif(1, 0, 1),
      both     = runif(1, 0, 4),
      1
    )
    curve[i] <- curve[i] * scale_factor
    curve
  }
  ensure_dir_exists <- function(path) {
    d <- dirname(path)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # Seed & setup --------------------------------------------------------------
  set.seed(seed)
  environmental_factors <- seq(env_min, env_max, length.out = env_length)
  cov_matrix <- build_cov_matrix(VE, VS, rES)

  # Simulate ------------------------------------------------------------------
  individual_norms <- do.call(rbind, lapply(seq_len(num_genotypes), function(i) {
    # random effects for replicates
    re <- MASS::mvrnorm(individuals_per_genotype, mu = c(0, 0), Sigma = cov_matrix)
    # fixed (mother) curve
    fixed_effect <- generate_fixed_norm(
      intercept = runif(1, 5, 30),
      slope     = runif(1, -1, 1),
      env       = environmental_factors,
      outliers  = outliers
    )
    # replicates
    reps <- do.call(rbind, lapply(seq_len(individuals_per_genotype), function(j) {
      u0 <- re[j,1]; u1 <- re[j,2]
      dev <- u0 + u1 * environmental_factors
      df <- data.frame(
        Genotype      = i,
        Replicate     = j,
        Environment   = environmental_factors,
        Trait         = fixed_effect + dev,
        ReactionNorm  = "Linear",
        BaseShift     = fixed_effect[1] + u0,
        Slope         = (fixed_effect[2] - fixed_effect[1]) /
          (environmental_factors[2] - environmental_factors[1]) + u1,
        VarianceOffset= VE,
        VarianceSlope = VS,
        Covariance    = rES * sqrt(VE * VS)
      )
      df
    }))
    # mother curve row
    mother <- data.frame(
      Genotype      = i,
      Replicate     = 0,
      Environment   = environmental_factors,
      Trait         = fixed_effect,
      ReactionNorm  = "Linear",
      BaseShift     = fixed_effect[1],
      Slope         = (fixed_effect[2] - fixed_effect[1]) /
        (environmental_factors[2] - environmental_factors[1]),
      VarianceOffset= NA,
      VarianceSlope = NA,
      Covariance    = NA
    )
    rbind(mother, reps)
  }))

  # Plotting ------------------------------------------------------------------
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    p <- {
      sel <- sample(unique(individual_norms$Genotype), genotypes_per_plot)
      d0  <- subset(individual_norms, Genotype %in% sel & Replicate == 0)
      d1  <- subset(individual_norms, Genotype %in% sel & Replicate != 0)
      ggplot2::ggplot() +
        ggplot2::geom_line(
          data = d0,
          ggplot2::aes(Environment, Trait, color = factor(Genotype), group = Genotype),
          size = 1
        ) +
        ggplot2::geom_line(
          data = d1,
          ggplot2::aes(Environment, Trait, color = factor(Genotype),
                       group = interaction(Genotype, Replicate)),
          linetype = "dashed", size = 0.5
        ) +
        ggplot2::labs(
          title = "Reaction Norms with Individual Deviations",
          x = "Environment", y = "Trait Value", color = "Genotype"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right")
    }
    if (!is.null(plot_file)) {
      ensure_dir_exists(plot_file)
      ggplot2::ggsave(plot_file, plot = p, width = 8, height = 6)
    }
  }

  # CSV output ----------------------------------------------------------------
  if (!is.null(csv_file)) {
    ensure_dir_exists(csv_file)
    utils::write.csv(individual_norms, file = csv_file, row.names = FALSE)
  }

  invisible(individual_norms)
}


