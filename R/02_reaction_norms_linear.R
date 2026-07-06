# =============================================================================
# 02_reaction_norms_linear.R — generate linear reaction norms
# =============================================================================
# WHAT IT DOES: Builds linear (intercept + slope) reaction norms across environments
#   for each genotype, optionally driven by the simulated genetics from 01.
# REQUIRES:     Genetic parameters from 01_simulate_genetics.R when USE_GENETICS = TRUE
#               (typically sourced upstream by 03_plasticity_scores.R).
# PRODUCES:     In memory: the linear reaction-norm data consumed by 03.
# HOW TO RUN:   Sourced by 03_plasticity_scores.R; or, with 01 already run,
#               setwd("~/PlastQuant"); source(here::here("R", "02_reaction_norms_linear.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit here)
#   USE_GENETICS  logical, default FALSE   use simulated genetics (FALSE = random params)  [COMMON]
#   SEED          integer, default 42      RNG seed                                        [COMMON]
#   OUTPUT_BASE   path, default here::here("output","default")  where outputs go
# -----------------------------------------------------------------------------
# These parameters are supplied by the sourcing script (03_plasticity_scores.R) or your
# R session. They are intentionally NOT assigned here: the script branches on whether
# SEED / USE_GENETICS exist (USE_GENETICS falls back to FALSE further down; SEED, when
# defined, offsets the RNG seed), so defining them here would change behaviour.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

# Function to build covariance matrix
build_cov_matrix <- function(VE, VS, rES) {
  matrix(c(
    VE, rES * sqrt(VE * VS),
    rES * sqrt(VE * VS), VS
  ), nrow = 2)
}


generate_fixed_norm <- function(intercept, slope, env, outliers = "None") {
  curve <- intercept + slope * env
  curve[curve < 0] <- 0
  if (outliers == "None") {
    return(curve)
  } else if (outliers == "positive") {
    i <- runif(1, 1, length(env))
    print("original")
    print(curve[i])
    curve[i] <- curve[i] * runif(1, 1, 2)
    print("mod")
    print(curve[i])
  } else if (outliers == "negative") {
    i <- runif(1, 1, length(env))
    print("original")
    print(curve[i])
    curve[i] <- curve[i] * runif(1, 0, 1)
    print("mod")
    print(curve[i])
  } else if (outliers == "both") {
    i <- runif(1, 1, length(env))
    print("original")
    print(curve[i])
    curve[i] <- curve[i] * runif(1, 0, 4)
    print("mod")
    print(curve[i])
  }

  return(curve)
}


generate_reaction_norm_with_covariance <- function(intercept, slope, env, cov_matrix) {
  random_effects <- MASS::mvrnorm(3, mu = c(0, 0), Sigma = cov_matrix)
  deviations <- t(sapply(1:3, function(i) {
    u0 <- random_effects[i, 1]
    u1 <- random_effects[i, 2]
    intercept_indiv <- intercept + u0
    slope_indiv <- slope + u1
    curve <- intercept_indiv + slope_indiv * env
    return(c(intercept_indiv, slope_indiv, curve))
  }))
  return(deviations)
}

introduce_outliers <- function() {}

create_plot <- function(individual_data, genotypes_per_plot = 5, title = "Reaction Norms with Individual Deviations") {
  selected_genotypes <- sample(unique(individual_data$Genotype), genotypes_per_plot)
  
  plot_data <- individual_data[individual_data$Genotype %in% selected_genotypes, ]
  
  ggplot() +
    geom_line(
      data = plot_data[plot_data$Replicate == 0, ],
      aes(x = Environment, y = Trait, color = factor(Genotype), group = Genotype),
      size = 1
    ) +
    geom_line(
      data = plot_data[plot_data$Replicate != 0, ],
      aes(x = Environment, y = Trait, color = factor(Genotype), group = interaction(Genotype, Replicate)),
      linetype = "dashed", size = 0.5
    ) +
    labs(title = title, x = "Environment", y = "Trait Value", color = "Genotype") + # Add a legend title for Genotype
    theme_minimal() +
    theme(legend.position = "right") # Ensure the legend is displayed on the right
}

create_individual_norm <- function(i, VE_vec, VS_vec, rES, num_individuals_per_genotype,
                                   fixed_intercepts, fixed_slopes, environmental_factors) {
  # Use genotype-specific VE and VS
  VE_i <- VE_vec[i]
  VS_i <- VS_vec[i]
  cov_matrix <- build_cov_matrix(VE_i, VS_i, rES)
  
  if (any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
    print(paste("Genotype:", i))
    print(paste("VE:", VE_i, "VS:", VS_i))
    print(cov_matrix)
  }
  
  # Sample random effects for the replicates for this genotype
  random_effects <- MASS::mvrnorm(num_individuals_per_genotype, mu = c(0, 0), Sigma = cov_matrix)
  
  # Compute the fixed (mother) curve with outlier modifications
  fixed_effect <- generate_fixed_norm(fixed_intercepts[i], fixed_slopes[i],
                                      environmental_factors,
                                      outliers = "none"
  )
  replicate_data <- do.call(rbind, lapply(1:num_individuals_per_genotype, function(j) {
    u0 <- random_effects[j, 1]
    u1 <- random_effects[j, 2]
  
    # Compute the additive deviation for each environment
    deviation <- u0 + u1 * environmental_factors
    
    # The replicate curve is the mother curve plus the individual deviation
    replicate_curve <- fixed_effect + deviation
    
    data.frame(
      Genotype = i,
      Replicate = j,
      Environment = environmental_factors,
      Trait = replicate_curve,
      ReactionNorm = "Linear",
      BaseShift = fixed_intercepts[i] + u0,
      Slope = fixed_slopes[i] + u1,
      VarianceOffset = VE_i,
      VarianceSlope = VS_i,
      Covariance = rES * sqrt(VE_i * VS_i)
    )
  }))

# Fixed data row for the mother curve
fixed_data <- data.frame(
  Genotype = i,
  Replicate = 0,
  Environment = environmental_factors,
  Trait = fixed_effect,
  ReactionNorm = "Linear",
  BaseShift = fixed_intercepts[i],
  Slope = fixed_slopes[i],
  VarianceOffset = NA,
  VarianceSlope = NA,
  Covariance = NA
)

# Combine the mother curve with the replicate data
return(rbind(fixed_data, replicate_data))
}

# -----------------------------------------------------------------------------
# EXECUTION
# -----------------------------------------------------------------------------
library(MASS)
library(ggplot2)

create_and_plot_linear_norms <- function(num_genotypes) {
  
  if (exists("SEED")) {
    seed_add_on = SEED;
  } else {
    seed_add_on = 0;
  }
  
  set.seed(42 + seed_add_on)
  environmental_factors <- seq(1, 10, length.out = 50) # 50 environments
  # Determine number of genotypes per form dynamically
  num_genotypes_per_form <- num_genotypes / 4
  
  # Use num_genotypes_per_form for generation
  num_genotypes_local <- num_genotypes_per_form
  num_individuals_per_genotype <- 3 # Number of individuals per genotype
  
  # -----------------------------------------------------------------------------
  # LOAD GENETIC PARAMETERS (OPTIONAL)
  # -----------------------------------------------------------------------------
  # Default to FALSE if not defined globally
  if (!exists("USE_GENETICS")) USE_GENETICS <- FALSE
  
  params_file <- paste(OUTPUT_BASE, "/synthetic_data/genetics/genotypic_parameters.csv",sep="")
  use_genetics_loaded <- FALSE
  
  # Initialize vectors
  fixed_intercepts <- numeric(num_genotypes_local)
  fixed_slopes <- numeric(num_genotypes_local)
  VE_vec <- rep(0.05, num_genotypes_local) # Default VE
  VS_vec <- rep(0.05, num_genotypes_local) # Default VS
  
  if (USE_GENETICS && file.exists(params_file)) {
    message("Loading parameters from genetic simulation: ", params_file)
    params_df <- read.csv(params_file)
  
    # Filter for Linear form (first quarter)
    if (nrow(params_df) >= num_genotypes_local) {
      params_subset <- params_df[1:num_genotypes_local, ]
      fixed_intercepts <- params_subset$BaseShift
      fixed_slopes <- params_subset$Slope
  
      if ("VE" %in% colnames(params_subset)) VE_vec <- params_subset$VE
      if ("VS" %in% colnames(params_subset)) VS_vec <- params_subset$VS
  
      use_genetics_loaded <- TRUE
    } else {
      warning("Not enough genotypes in parameter file. Falling back to random generation.")
    }
  }
  
  if (!use_genetics_loaded) {
    if (USE_GENETICS) message("Genetic parameter file not found or invalid. Using random generation.")
    fixed_intercepts <- runif(num_genotypes_local, 5, 30)
    fixed_slopes <- runif(num_genotypes_local, -1, 1)
  }
  
  # Covariance matrix setup (Global rES)
  rES <- 0.5
  
  # Generate individual (replicate) norms for each genotype
  
  individual_norms <- do.call(rbind, lapply(1:num_genotypes_local, create_individual_norm, 
                                            VE_vec, VS_vec, rES, num_individuals_per_genotype,
                                          fixed_intercepts, fixed_slopes, environmental_factors))
  
  
  p <- create_plot(individual_norms, genotypes_per_plot = 3, title = "Random Sample of 3 Genotypes")
  pdf(paste(OUTPUT_BASE, "/plots/linear_reaction_norms2.pdf",sep=""), width = 8, height = 6)
  print(p)
  dev.off()
  output_file <- paste(OUTPUT_BASE, "/synthetic_data/fixed_full/linear_reaction_norms_data.csv",sep="")
  
  ensure_dir_exists <- function(file_path) {
    dir_path <- dirname(file_path)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  ensure_dir_exists(output_file)
  
  
  # write out without row names
  write.csv(individual_norms, file = output_file, row.names = FALSE)
  return(individual_norms)
}
