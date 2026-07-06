# =============================================================================
# 02_reaction_norms_nonlinear.R — generate nonlinear reaction norms
# =============================================================================
# WHAT IT DOES: Builds nonlinear reaction norms (Gaussian, sinusoidal, and wave
#   shapes) across environments for each genotype, optionally driven by the
#   simulated genetics from 01.
# REQUIRES:     Genetic parameters from 01_simulate_genetics.R when USE_GENETICS = TRUE
#               (typically sourced upstream by 03_plasticity_scores.R).
# PRODUCES:     In memory: the nonlinear reaction-norm data consumed by 03.
# HOW TO RUN:   Sourced by 03_plasticity_scores.R; or, with 01 already run,
#               setwd("~/PlastQuant"); source(here::here("R", "02_reaction_norms_nonlinear.R"))
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

library(MASS) # For multivariate normal distributions
library(ggplot2) # For plotting
library(dplyr) # For data manipulation

# NON-LINEAR REACTION NORM FUNCTIONS
generate_gaussian_reaction_norm <- function(env, base_shift, slope, amplitude = 4, width = 0.2, center = 5) {
  base_shift + slope * amplitude * exp(-width * (env - center)^2)
}

generate_sinusoidal_reaction_norm <- function(env, base_shift, slope, amplitude = 4, frequency = 0.5, phase = 0) {
  base_shift + slope * amplitude * sin(frequency * env + phase)
}

generate_wave_reaction_norm <- function(env, base_shift, slope, amplitude = 7, frequency = 1.2) {
  base_shift + slope * amplitude * sin(frequency * (env - 5)) * exp(-0.1 * (env - 5)^2)
}

# Function to build covariance matrix
build_cov_matrix <- function(var_offset, var_slope, correlation) {
  cov <- correlation * sqrt(var_offset * var_slope)
  matrix(c(var_offset, cov, cov, var_slope), nrow = 2)
}

# New function to generate a fixed non-linear norm with outlier modifications
generate_fixed_nonlin_norm <- function(env, base_shift, slope, norm_function, outliers = "None") {
  base_curve <- norm_function(env, base_shift, slope)
  
  if (outliers == "None") {
    return(base_curve)
  } else if (outliers == "positive") {
    i <- sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] <- base_curve[i] * runif(1, 1, 2)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  } else if (outliers == "negative") {
    i <- sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] <- base_curve[i] * runif(1, 0, 1)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  } else if (outliers == "both") {
    i <- sample(1:length(env), 1)
    cat("Original value at index", i, ":", base_curve[i], "\n")
    base_curve[i] <- base_curve[i] * runif(1, 0, 4)
    cat("Modified value at index", i, ":", base_curve[i], "\n")
  }
  
  return(base_curve)
}

# Modified function to generate data for a single genotype that inherits the mother's outlier
generate_genotype_data <- function(base_shift, slope, env, norm_function, num_replicates,
                                   var_offset, var_slope, correlation, genotype_id, norm_type,
                                   outliers = "None") {
  cov_matrix <- build_cov_matrix(var_offset, var_slope, correlation)
  
  # Compute the mother (fixed) curve with outlier modifications
  mother_curve <- generate_fixed_nonlin_norm(env, base_shift, slope, norm_function, outliers)
  # Also compute the base curve (without outliers) for calculating deviations
  base_curve <- norm_function(env, base_shift, slope)
  
  mother_data <- data.frame(
    Genotype = genotype_id,
    Replicate = 0,
    Environment = env,
    Trait = mother_curve,
    ReactionNorm = norm_type,
    BaseShift = base_shift,
    Slope = slope,
    VarianceOffset = NA,
    VarianceSlope = NA,
    Covariance = NA
  )
  
  # Generate replicates that inherit the mother's outlier pattern
  replicate_data <- do.call(rbind, lapply(1:num_replicates, function(rep) {
    perturbations <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = cov_matrix)
    offset_perturbation <- perturbations[1]
    slope_perturbation <- perturbations[2]
    
    # Compute the perturbed curve without outlier modifications
    perturbed_curve <- norm_function(env, base_shift + offset_perturbation, slope + slope_perturbation)
    # Calculate the deviation relative to the base (non-outlier) curve
    deviation <- perturbed_curve - base_curve
    # Add the deviation to the mother curve so that the outlier is inherited
    replicate_curve <- mother_curve + deviation
    
    data.frame(
      Genotype = genotype_id,
      Replicate = rep,
      Environment = env,
      Trait = replicate_curve,
      ReactionNorm = norm_type,
      BaseShift = base_shift + offset_perturbation,
      Slope = slope + slope_perturbation,
      VarianceOffset = var_offset,
      VarianceSlope = var_slope,
      Covariance = correlation * sqrt(var_offset * var_slope)
    )
  }))
  
  return(rbind(mother_data, replicate_data))
}

# Modified function to generate data for all genotypes for a specific norm type
generate_reaction_norm_data <- function(norm_function, norm_type, base_shifts, slopes,
                                        env, num_replicates, var_offset, var_slope, correlation,
                                        amplitudes = NULL, widths = NULL, centers = NULL, frequencies = NULL, phases = NULL,
                                        outliers = "None", id_offset = 0) {
  # Generate IDs based on offset
  genotype_ids <- (id_offset + 1):(id_offset + length(base_shifts))
  
  do.call(rbind, lapply(seq_along(genotype_ids), function(i) {
    genotype_id <- genotype_ids[i]
    base_shift <- base_shifts[i]
    slope <- slopes[i]
    
    # Handle vector inputs for variance parameters
    vo <- if (length(var_offset) > 1) var_offset[i] else var_offset
    vs <- if (length(var_slope) > 1) var_slope[i] else var_slope
    
    # Handle shape-specific parameters (if provided)
    # Note: The generator functions (generate_gaussian_reaction_norm etc.) need to be updated
    # or we need to pass these as arguments.
    # Currently, the generator functions have defaults. We need to pass them explicitly.
    # But `generate_genotype_data` calls `norm_function(env, base_shift, slope)`.
    # We need to modify `generate_genotype_data` or use a wrapper.
    # Actually, `generate_genotype_data` is defined in this file (lines 66-107).
    # It calls `norm_function(env, base_shift, slope)`.
    # We need to pass the extra args to `norm_function`.
    
    # Let's construct the arguments for the norm function dynamically
    args_list <- list(env = env, base_shift = base_shift, slope = slope)
    if (!is.null(amplitudes)) args_list$amplitude <- amplitudes[i]
    if (!is.null(widths)) args_list$width <- widths[i]
    if (!is.null(centers)) args_list$center <- centers[i]
    if (!is.null(frequencies)) args_list$frequency <- frequencies[i]
    if (!is.null(phases)) args_list$phase <- phases[i]
    
    # We need to modify generate_genotype_data to accept these extra args or pass the function closure.
    # Simpler approach: Create a closure for the norm function with fixed parameters for this genotype.
    specific_norm_func <- function(e, b, s) {
      # We ignore b and s here because we capture the specific ones,
      # BUT generate_genotype_data passes them.
      # Better: do.call the original norm_function with all args.
      do.call(norm_function, c(list(env = e, base_shift = b, slope = s), args_list[-(1:3)]))
    }
    
    generate_genotype_data(
      base_shift, slope, env, specific_norm_func, num_replicates,
      vo, vs, correlation, genotype_id, norm_type, outliers
    )
  }))
}

# Plotting functions (unchanged)
create_plot <- function(data, genotypes_per_plot = 3, show_replicates = TRUE) {
  if (exists("SEED")) {
    seed_add_on = SEED;
  } else {
    seed_add_on = 0;
  }
  set.seed(44 + seed_add_on)
  selected_genotypes <- sample(unique(data$Genotype), genotypes_per_plot, replace = TRUE)
  plot_data <- data[data$Genotype %in% selected_genotypes, ]
  
  ggplot(plot_data, aes(x = Environment, y = Trait, color = factor(Genotype))) +
    geom_line(aes(group = interaction(Genotype, Replicate), linetype = as.factor(Replicate))) +
    facet_wrap(~ReactionNorm) +
    labs(
      title = "Reaction Norms",
      x = "Environment",
      y = "Trait Value",
      color = "Genotype",
      linetype = "Replicate"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}

create_plot2 <- function(gaussian_data, sinusoidal_data, wave_data, show_replicates = TRUE) {
  if (exists("SEED")) {
    seed_add_on = SEED;
  } else {
    seed_add_on = 0;
  }
  set.seed(44 + seed_add_on)
  
  gaussian_sample <- gaussian_data %>%
    filter(Genotype == sample(unique(Genotype), 1))
  
  sinusoidal_sample <- sinusoidal_data %>%
    filter(Genotype == sample(unique(Genotype), 1))
  
  wave_sample <- wave_data %>%
    filter(Genotype == sample(unique(Genotype), 1))
  
  sampled_data <- rbind(gaussian_sample, sinusoidal_sample, wave_sample)
  
  fixed_data <- sampled_data[sampled_data$Replicate == 0, ]
  replicate_data <- sampled_data[sampled_data$Replicate != 0, ]
  
  p <- ggplot(fixed_data, aes(x = Environment, y = Trait, color = ReactionNorm, group = interaction(ReactionNorm, Genotype))) +
    geom_line(size = 1, linetype = "solid") +
    labs(
      title = "Randomly Sampled Reaction Norms (One Per Type)",
      x = "Environment",
      y = "Trait Value",
      color = "Reaction Norm Type"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  if (show_replicates) {
    p <- p + geom_line(
      data = replicate_data,
      aes(group = interaction(ReactionNorm, Genotype, Replicate)),
      linetype = "dashed", size = 0.5
    )
  }
  
  return(p)
}

create_and_plot_nonlinear_norms <- function(num_genotypes = 200) {

  if (exists("SEED")) {
    seed_add_on = SEED;
  } else {
    seed_add_on = 0;
  }
  set.seed(12 + seed_add_on)
  environmental_factors <- seq(1, 10, length.out = 50) # 50 environments
  # Determine number of genotypes per form dynamically
  if (!exists("num_genotypes")) num_genotypes <- 200 # Default if not set
  num_genotypes_per_form <- num_genotypes / 4
  
  num_replicates <- 3 # Number of replicates for each genotype
  
  
  # -----------------------------------------------------------------------------
  # EXECUTION
  # -----------------------------------------------------------------------------
  # -----------------------------------------------------------------------------
  # LOAD GENETIC PARAMETERS (OPTIONAL)
  # -----------------------------------------------------------------------------
  # Default to FALSE if not defined globally
  if (!exists("USE_GENETICS")) USE_GENETICS <- FALSE
  
  params_file <- paste(OUTPUT_BASE, "/synthetic_data/genetics/genotypic_parameters.csv",sep="")
  use_genetics_loaded <- FALSE
  
  # Initialize parameter lists with defaults
  # (Using vectors of length num_genotypes_per_form)
  params_gaussian   <- list(shifts = runif(num_genotypes_per_form, 10, 30), slopes = runif(num_genotypes_per_form, 0.5, 1.5), ve = 0.05, vs = 0.01, amp = rep(4,   num_genotypes_per_form), width  = rep(0.2, num_genotypes_per_form), center = rep(5,   num_genotypes_per_form))
  params_sinusoidal <- list(shifts = runif(num_genotypes_per_form, 10, 30), slopes = runif(num_genotypes_per_form, 0.5, 1.5), ve = 0.05, vs = 0.01, amp = rep(4,   num_genotypes_per_form), freq   = rep(0.5, num_genotypes_per_form), phase  = rep(0,   num_genotypes_per_form))
  params_wave       <- list(shifts = runif(num_genotypes_per_form, 10, 30), slopes = runif(num_genotypes_per_form, 0.5, 1.5), ve = 0.05, vs = 0.01, amp = rep(7,   num_genotypes_per_form), freq   = rep(1.2, num_genotypes_per_form))
  
  
  if (USE_GENETICS && file.exists(params_file)) {
    message("Loading parameters from genetic simulation: ", params_file)
    params_df <- read.csv(params_file)
  
    if (nrow(params_df) >= num_genotypes) {
      # Helper to extract
      extract_params <- function(df_subset) {
        p <- list(
          shifts = df_subset$BaseShift,
          slopes = df_subset$Slope,
          ve = if ("VE" %in% colnames(df_subset)) df_subset$VE else 0.5,
          vs = if ("VS" %in% colnames(df_subset)) df_subset$VS else 0.1,
          amp = if ("Amplitude" %in% colnames(df_subset)) df_subset$Amplitude else 4,
          width = if ("Width" %in% colnames(df_subset)) df_subset$Width else 0.2,
          center = if ("Center" %in% colnames(df_subset)) df_subset$Center else 5,
          freq = if ("Frequency" %in% colnames(df_subset)) df_subset$Frequency else 0.5,
          phase = if ("Phase" %in% colnames(df_subset)) df_subset$Phase else 0
        )
        return(p)
      }
  
      # Dynamic Ranges
      n <- num_genotypes_per_form
      params_gaussian <- extract_params(params_df[(n + 1):(2 * n), ])
      params_sinusoidal <- extract_params(params_df[(2 * n + 1):(3 * n), ])
      params_wave <- extract_params(params_df[(3 * n + 1):(4 * n), ])
      print("Default Gaussian Parameters:")
      print(params_gaussian)
      use_genetics_loaded <- TRUE
    } else {
      warning(paste("Not enough genotypes in parameter file (need", num_genotypes, "). Falling back to random generation."))
    }
  }
  
  # Generate Data
  gaussian_data <- generate_reaction_norm_data(
    generate_gaussian_reaction_norm, "Gaussian", params_gaussian$shifts, params_gaussian$slopes,
    environmental_factors, num_replicates,
    var_offset = params_gaussian$ve, var_slope = params_gaussian$vs, correlation = -0.5,
    amplitudes = params_gaussian$amp, widths = params_gaussian$width, centers = params_gaussian$center,
    outliers = "none", id_offset = num_genotypes_per_form
  )
  sinusoidal_data <- generate_reaction_norm_data(
    generate_sinusoidal_reaction_norm, "Sinusoidal", params_sinusoidal$shifts, params_sinusoidal$slopes,
    environmental_factors, num_replicates,
    var_offset = params_sinusoidal$ve, var_slope = params_sinusoidal$vs, correlation = -0.3,
    amplitudes = params_sinusoidal$amp, frequencies = params_sinusoidal$freq, phases = params_sinusoidal$phase,
    outliers = "none", id_offset = 2 * num_genotypes_per_form
  )
  wave_data <- generate_reaction_norm_data(
    generate_wave_reaction_norm, "Wave", params_wave$shifts, params_wave$slopes,
    environmental_factors, num_replicates,
    var_offset = params_wave$ve, var_slope = params_wave$vs, correlation = -0.7,
    amplitudes = params_wave$amp, frequencies = params_wave$freq,
    outliers = "none", id_offset = 3 * num_genotypes_per_form
  )
  
  all_data <- rbind(gaussian_data, sinusoidal_data, wave_data)
  
  
  
  pdf(paste(OUTPUT_BASE, "/plots/reaction_norms_plots_combined.pdf",sep=""), width = 12, height = 8)
  print(create_plot(gaussian_data)) # Plot for Gaussian norms
  print(create_plot(sinusoidal_data)) # Plot for Sinusoidal norms
  print(create_plot(wave_data)) # Plot for Wave norms
  print(create_plot2(gaussian_data, sinusoidal_data, wave_data)) # Combined plot with 3 randomly sampled norms
  dev.off()
  
  ensure_dir_exists <- function(file_path) {
    dir_path <- dirname(file_path)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  ensure_dir_exists(paste(OUTPUT_BASE, "/synthetic_data/fixed_full",sep=""))
  write.csv(gaussian_data, paste(OUTPUT_BASE, "/synthetic_data/fixed_full/gaussian_data.csv",sep=""), row.names = FALSE)
  write.csv(sinusoidal_data, paste(OUTPUT_BASE, "/synthetic_data/fixed_full/sinusoidal_data.csv",sep=""), row.names = FALSE)
  write.csv(wave_data, paste(OUTPUT_BASE, "/synthetic_data/fixed_full/wave_data.csv",sep=""), row.names = FALSE)
  write.csv(all_data, paste(OUTPUT_BASE, "/synthetic_data/fixed_full/all_combined_data.csv",sep=""), row.names = FALSE)

  return(list(gaussian_data = gaussian_data, sinusoidal_data = sinusoidal_data, wave_data = wave_data))
}
