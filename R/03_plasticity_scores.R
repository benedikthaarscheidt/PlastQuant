# =============================================================================
# 03_plasticity_scores.R — compute every plasticity index (pipeline hub)
# =============================================================================
# WHAT IT DOES: Computes every plasticity index (via the ppindices package) for each
#   simulated genotype's reaction norm across the four forms (linear, gaussian,
#   sinusoidal, wave), assembling the `scores_list` consumed by GWAS and analysis.
#   This is the hub: it auto-sources 01 and 02, so running it reruns the upstream
#   simulation.
# REQUIRES:     Auto-sources 01_simulate_genetics.R, 02_reaction_norms_linear.R,
#               02_reaction_norms_nonlinear.R; library(ppindices).
#               When USE_GENETICS = TRUE (the default) the caller MUST set
#               HERITABILITY_5TH, CAUSAL_SNP_NUM and OUTPUT_BASE first (a scenario_*.R
#               driver, or assignments in your R session) — the script stops otherwise.
# PRODUCES:     `scores_list` in memory; per-score CSVs under OUTPUT_BASE when SAVE=TRUE.
# HOW TO RUN:   setwd("~/PlastQuant")
#               HERITABILITY_5TH <- 1; CAUSAL_SNP_NUM <- 5
#               OUTPUT_BASE <- here::here("output","my_run")
#               source(here::here("R", "03_plasticity_scores.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these in a scenario_*.R driver or your R session before sourcing.
#   The defaults below are applied in-script (guarded) at the noted line numbers.
#   USE_GENETICS          logical, default TRUE  use simulated genetics (line 36)        [COMMON]
#   HERITABILITY_5TH      numeric, REQUIRED when USE_GENETICS=TRUE  H^2 of the 5th trait [COMMON]
#   CAUSAL_SNP_NUM        integer, REQUIRED when USE_GENETICS=TRUE  causal SNPs / param  [COMMON]
#   OUTPUT_BASE           path,    REQUIRED when USE_GENETICS=TRUE  output location      [COMMON]
#   NUM_GENOTYPES         integer, default 800   number of genotypes (line 100)
#   N_MEAS_POINTS         integer, default 6     measurement points per reaction norm (line 161)
#   KEEP_REPLICATES       logical, default (see lines 97/179)  keep replicate rows
#   STRUCTURED_POPULATION logical, default FALSE nested population structure (line 98)
#   GENETIC_VARIANCES     logical/vector         per-parameter genetic variances (from 01)
#   POLY_MODE/STRENGTH/COUNTS  causal-effect controls passed through to 01
#   GWAS                  logical                run GWAS downstream (read by 04/05/06)
#   SAVE                  logical                write per-score CSVs under OUTPUT_BASE
# =============================================================================
# Original notes:
# Script:    03_plasticity_scores.R
# Purpose:   Load all plasticity‐score functions, set sampling parameters (range,
#            interval, indices, covariate, type labels, etc.), preprocess reaction‐
#            norm data into per‐genotype trait vectors, and compute every plasticity
#            score across the four genotype forms (linear, gaussian, sinusoidal, wave).
#
# Steps:
#   1. Source the five score‐definition files (1_2.R, 2.R, 3.R, 4.R, 5.R) plus
#      the linear and nonlinear norms generators.
#   2. Define default globals (indices, interval, env, Covariate, type_labels) and
#      only override them if already defined in the calling environment (which happens when the script is launched from the resolution.R file).
#   3. Ensure the raw reaction‐norm data frames (gaussian_data, sinusoidal_data,
#      wave_data, individual_norms) are present, assigning defaults if needed.
#   4. Combine the four raw data sets into one data frame for any combined‐data ops.
#   5. Preprocess each form’s data via `preprocess_data()` to:
#        • average replicates per genotype × environment,
#        • subsample at `indices`,
#        • return a named list of genotype‐vectors.
#   6. For each of the plasticity functions :
#        • Call them over each genotype list via `call_function()`,
#        • Post‐process into a unified data frame with `Type` labels.
#   7. Build four reaction‐norm matrices (`linear`, `gaussian`, `sinusoidal`, `wave`)
#      by row‐binding the per‐genotype vectors.
#   8. Make “long” trait + environment data frames for plotting or summary stats.
#
# Outputs (in memory):
#   • `processed_data_list`: list of four preprocessed genotype‐lists
#   • Score tables
#   • Reaction‐norm matrices: `linear_matrix`, `gaussian_matrix`, `sinusoidal_matrix`,
#     `wave_matrix`
#
# Globals used/produced:
#   indices, interval, env, Covariate, type_labels,
#   gaussian_data, sinusoidal_data, wave_data, linear_data

if (!exists("USE_GENETICS")) USE_GENETICS <- TRUE

if (USE_GENETICS) {
  if (!exists("HERITABILITY_5TH") || !exists("CAUSAL_SNP_NUM") || !exists("OUTPUT_BASE")) {
    stop("HERITABILITY_5TH, CAUSAL_SNP_NUM and OUTPUT_BASE need to be set in the sourcing script")
  }
} else {
  if (!exists("OUTPUT_BASE")) OUTPUT_BASE <- here::here("output", "no_genetics_output")
  if (!exists("HERITABILITY_5TH")) HERITABILITY_5TH <- 0
  if (!exists("CAUSAL_SNP_NUM")) CAUSAL_SNP_NUM <- 0
}

if (!dir.exists(paste(OUTPUT_BASE,"/plots",sep=""))) dir.create(paste(OUTPUT_BASE,"/plots",sep=""), recursive = TRUE)

library(ppindices)
# coefficient-of variation total (calculate_CVt) - tested,
# slope of norm reaction (calculate_reaction_norm_slope) - tested,
# slope of plastic response (D) (calculate_D_slope)- tested,
# response coefficient (RC) (calculate_RC)- tested,
# Standard deviation of means (CVm) (calculate_CVm)- tested,
# Standard deviation of medians (CVmd)(calculate_CVmd)- tested,
# Grand plasticity (calculate_GPi)- tested,
# Phenotypic Plasticity Index (calculate_PPF)- tested,
# Phenotypic Plasticity Index (calculate_Phenotypic_Plasticity_Index)- tested,
# PImd (calculate_PImd)- tested,
# PILSM (calculate_PILSM)- tested,
# RTR (calculate_RTR)- tested,
# PIR (calculate_PIR) - tested
# RDPI	(rdpi_calculation) - tested,
# RDPIs (rdpi_mean_calculation) - tested,
# ESPI (calculate_ESPI) - tested,
# ESPIid (espiid_calculation) - tested,
# evwpi_calculation (idea from Benedikt)
# Phenotypic Stability Index (calculate_PSI),
# Relative Plasticity Index (calculate_RPI) - tested,
# Plasticity Quotient (calculate_PQ) - tested,
# Phenotypic Range (PR) (calculate_PR) - tested,
# Norm of reaction width (calculate_NRW) - tested,
# Environment-Specific Plasticity (ESP) (calculate_ESP) - tested,
# Calculate Plasticity Differential (PD) (calculate_PD) - tested,
# Calculate Fitness Plasticity Index (FPI) (calculate_FPI) - tested,
# Calculate Transplant Plasticity Score (TPS)(calculate_TPS) - tested,
# Calculate Developmental Plasticity Index (DPI)(calculate_DPI) - tested,
# Calculate Coefficient of Environmental Variation (CEV)(calculate_CEV) - tested,
# Calculate Plasticity Response Index (PRI)(calculate_PRI) - tested,
# Calculate Phenotypic Flexibility Index (PFI)(calculate_PFI) - tested,
# Calculate Standardized Plasticity Index (SPI)(calculate_SPI) - tested,
# Calculate Absolute Plasticity Coefficient (APC)(calculate_APC) - tested,
# Calculate Stability Index (SI)(calculate_SI) - tested,
# Calculate Relative Stability Index (RSI)(calculate_RSI),
# Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
# Calculate Environmental Variance Sensitivity (EVS)(calculate_EVS) - tested,
# Calculate Multivariate Plasticity Index (MVPi)(calculate_MVPi) NEEDS TO BE TESTED WITH THE REQUESTED DATASET FROM PROF. BARBOSA,
# Calculate Standardized Plasticity Metric (SPM)(calculate_SPM) - tested,
# Calculate SSpop/SStotal Plasticity Ratio(calculate_Plasticity_Ratio) - tested


# -----------------------------------------------------------------------------
# ENABLE GENETICS AND GENERATE DATA
# -----------------------------------------------------------------------------
# USE_GENETICS is already set above (default TRUE, or FALSE if set before sourcing)
if (!exists("KEEP_REPLICATES")) KEEP_REPLICATES <- FALSE
if (!exists("STRUCTURED_POPULATION")) STRUCTURED_POPULATION <- FALSE
if (!exists("GENETIC_VARIANCES")) GENETIC_VARIANCES <- FALSE
if (!exists("NUM_GENOTYPES")) NUM_GENOTYPES <- 800
if (!exists("NUM_SNPs")) NUM_SNPs <- 200000

# --- Simulation Configuration ---

if (USE_GENETICS) {
  # Total number of genotypes (must be divisible by 4)
  POLY_COUNTS <- list(
    Slope = CAUSAL_SNP_NUM,
    BaseShift = CAUSAL_SNP_NUM,
    Amplitude = CAUSAL_SNP_NUM,
    Width = CAUSAL_SNP_NUM,
    Center = CAUSAL_SNP_NUM,
    Frequency = CAUSAL_SNP_NUM,
    Phase = CAUSAL_SNP_NUM,
    VE = 0,
    VS = 0
  )

  # Scale 1-5: How strongly genetics affects the parameter (Scales the Genetic SD)
  # 5 = Full SD (High Heritability), 1 = 20% SD (Low Heritability)
  POLY_STRENGTH <- list(
    Slope = HERITABILITY_5TH,
    BaseShift = HERITABILITY_5TH,
    Amplitude = HERITABILITY_5TH,
    Width = HERITABILITY_5TH,
    Center = HERITABILITY_5TH,
    Frequency = HERITABILITY_5TH,
    Phase = HERITABILITY_5TH,
    VE = 0,
    VS = 0
  )

  # Effect Size Distribution:
  # "uniform"     = All SNPs have similar effect sizes (Box-like)
  # "exponential" = Few strong SNPs, many weak ones (L-shaped, more realistic)
  POLY_MODE <- "uniform"

  message(paste("Step 1: Simulating Genetics (", NUM_GENOTYPES, " Genotypes, Polygenic Model)..."))
  source(here::here("R", "01_simulate_genetics.R"))
  sim_gen_res = simulate_genetics()
  snp_matrix = sim_gen_res$snp_matrix
  params_df = sim_gen_res$params_df
  truth_df = sim_gen_res$truth_df
  
  message("Step 2: Generating Reaction Norms from Genetic Parameters...")
} else {
  message("Step 1: Skipping genetics (USE_GENETICS = FALSE). Reaction norms will use random parameters.")
  message("Step 2: Generating Reaction Norms with random parameters (no genetic influence)...")
}

source(here::here("R", "02_reaction_norms_linear.R"))
linear_data = create_and_plot_linear_norms(NUM_GENOTYPES)
source(here::here("R", "02_reaction_norms_nonlinear.R"))
nonlinear_data = create_and_plot_nonlinear_norms(NUM_GENOTYPES)
gaussian_data = nonlinear_data$gaussian_data
sinusoidal_data = nonlinear_data$sinusoidal_data
wave_data = nonlinear_data$wave_data

## ──────────────────────────────────────────────────────────────────────────────
##  1) Define defaults for any globals you need
if (!exists("N_MEAS_POINTS")) N_MEAS_POINTS <- 6
default_indices <- round(seq(1, 50, length.out = N_MEAS_POINTS)) # equidistant indices from 1 to 50
default_env <- seq(1, 10, length.out = length(default_indices))
default_Covariate <- rep(1, length(default_indices))
default_type_labels <- c("gaussian", "sinusoidal", "wave", "linear")


## ──────────────────────────────────────────────────────────────────────────────å
##  2) Only assign them if they don't already exist
## ──────────────────────────────────────────────────────────────────────────────
if (!exists("indices", envir = .GlobalEnv)) indices <- default_indices
if (!exists("env", envir = .GlobalEnv)) env <- default_env
if (!exists("Covariate", envir = .GlobalEnv)) Covariate <- default_Covariate
if (!exists("type_labels", envir = .GlobalEnv)) type_labels <- default_type_labels

message("Using indices: ", paste(indices, collapse = ", "))

# Default KEEP_REPLICATES to TRUE if not defined (for standalone execution)
if (!exists("KEEP_REPLICATES", envir = .GlobalEnv)) KEEP_REPLICATES <- TRUE

combined_df <- rbind(linear_data, gaussian_data, sinusoidal_data, wave_data)

# -----------------------------------------------------------------------------
# PLOTTING CHECK (Visual Verification)there
# -----------------------------------------------------------------------------


library(ggplot2)
library(gridExtra)

plot_genotypes <- function(df, form_name, n_plot = 3) {
  # Get unique genotypes
  genos <- unique(df$Genotype)
  # Pick random subset
  selected_genos <- sample(genos, min(length(genos), n_plot))

  subset_df <- df[df$Genotype %in% selected_genos, ]

  p <- ggplot(subset_df, aes(x = Environment, y = Trait, group = interaction(Genotype, Replicate), color = as.factor(Genotype))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Genotype, scales = "free_y") +
    labs(
      title = paste(form_name, "- Sample Genotypes"),
      y = "Phenotype", x = "Environment", color = "Genotype"
    ) +
    theme_minimal()
  return(p)
}

pdf(paste(OUTPUT_BASE,"/plots/reaction_norms_check.pdf",sep=""), width = 10, height = 8)
print(plot_genotypes(linear_data, "Linear"))
print(plot_genotypes(gaussian_data, "Gaussian"))
print(plot_genotypes(sinusoidal_data, "Sinusoidal"))
print(plot_genotypes(wave_data, "Wave"))
dev.off()

message(paste("Reaction norm plots saved to: ", OUTPUT_BASE,"/plots/reaction_norms_check.pdf",sep=""))


############################################################

#' Preprocess Dataset for Modular Index Functions
#'
#' This function preprocesses a dataset by averaging trait values over replicates
#' for each genotype and environmental condition. It outputs a list of trait vectors,
#' one for each genotype, with the length corresponding to the number of environmental conditions.
#'
#' @param data A data frame containing the dataset.
#' @param trait_col The column name or index containing the trait values.
#' @param genotype_col The column name or index containing the genotype information.
#' @param replicate_col (Optional) The column name or index containing replicate information.
#'   If provided, the function averages trait values across replicates for each genotype-environment combination.
#' @return A named list where each element corresponds to a genotype and contains:
#'   - A numeric vector of averaged trait values for that genotype.
#'
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   Genotype = rep(c("G1", "G2"), each = 10),
#'   Environment = rep(1:5, times = 4),
#'   Replicate = rep(1:2, each = 5, times = 2),
#'   Trait = c(
#'     10, 12, 15, 20, 22, 15, 17, 20, 18, 20,
#'     23, 25, 28, 30, 32, 25, 27, 30, 28, 30
#'   )
#' )
#'
#' # Preprocess the dataset
#' preprocessed <- preprocess_data(
#'   data,
#'   trait_col = "Trait", genotype_col = "Genotype",
#'   replicate_col = "Replicate"
#' )
#'
#' print(preprocessed)
#'
#' @export
preprocess_data <- function(data, trait_col, genotype_col, replicate_col = NULL, indices = NULL) {
  if (is.numeric(trait_col)) trait_col <- colnames(data)[trait_col]
  if (is.numeric(genotype_col)) genotype_col <- colnames(data)[genotype_col]
  if (!is.null(replicate_col) && is.numeric(replicate_col)) replicate_col <- colnames(data)[replicate_col]

  genotype_split <- split(data, data[[genotype_col]])
  result <- list()

  for (genotype in names(genotype_split)) {
    genotype_data <- genotype_split[[genotype]]

    if (is.null(replicate_col)) {
      result[[genotype]] <- genotype_data[[trait_col]][seq(1, length(genotype_data[[trait_col]]), 1)]
    } else {
      unique_replicates <- unique(genotype_data[[replicate_col]])

      temp <- do.call(
        cbind,
        lapply(unique_replicates, function(rep) {
          genotype_data[genotype_data[[replicate_col]] == rep, trait_col]
        })
      )

      averaged_vector <- rowMeans(temp, na.rm = TRUE)

      # If indices are provided externally, use them.
      # Otherwise, compute default indices over the entire vector.
      if (is.null(indices)) {
        indices <- unique(c(
          seq(1, length(averaged_vector), by = interval),
          length(averaged_vector)
        ))
      }

      if (exists("KEEP_REPLICATES") && KEEP_REPLICATES) {
        # Return each replicate as a separate entry in the list
        for (rep_idx in seq_along(unique_replicates)) {
          rep_val <- unique_replicates[rep_idx]
          # Extract the specific replicate vector
          rep_vector <- genotype_data[genotype_data[[replicate_col]] == rep_val, trait_col]
          # Apply indices subsampling
          result[[paste0(genotype, "_Rep_", rep_val)]] <- rep_vector[indices]
        }
      } else {
        # Default behavior: Average replicates
        result[[genotype]] <- averaged_vector[indices]
      }
    }
  }

  return(result)
}


call_function <- function(data_list, score, additional_1 = NULL, additional_2 = NULL) {
  results <- list() # Store results for each item

  for (i in names(data_list)) {
    # Construct arguments dynamically
    args <- list(data_list[[i]])
    if (!is.null(additional_1)) args <- c(args, list(additional_1))
    if (!is.null(additional_2)) args <- c(args, list(additional_2))

    # Call the score function
    result <- do.call(score, args)

    # Store the result in the list
    results[[i]] <- result
  }

  # Combine results into a data frame or keep as list
  if (all(sapply(results, is.list))) {
    # Extract all unique names from the output lists
    output_names <- unique(unlist(lapply(results, names)))

    # Create a data frame with one column per output name
    results_df <- do.call(rbind, lapply(results, function(res) {
      # Fill missing outputs with NA
      sapply(output_names, function(name) if (name %in% names(res)) res[[name]] else NA)
    }))

    rownames(results_df) <- names(results)
    return(as.data.frame(results_df))
  } else {
    # Ensure consistent return of numeric vector for scalar results
    return(unlist(results))
  }
}

post_process <- function(data_list, score_function, type_labels, ...) {
  result_list <- list()

  for (i in seq_along(data_list)) {
    # Calculate scores using the provided score function and additional arguments
    score_data <- call_function(data_list[[i]], score_function, ...)

    # Convert to dataframe if needed and add the Type label
    if (is.null(dim(score_data))) {
      # score_data is a vector or named list
      score_df <- data.frame(
        Score = as.numeric(score_data),
        Type = type_labels[i],
        stringsAsFactors = FALSE,
        row.names = names(score_data)
      )
    } else if (is.matrix(score_data)) {
      # score_data is a matrix
      score_df <- as.data.frame(score_data, stringsAsFactors = FALSE)
      score_df$Type <- type_labels[i]
      # Rename first column to "Score" if it exists
      if (ncol(score_df) > 1) {
        colnames(score_df)[1] <- "Score"
      } else {
        colnames(score_df) <- c("Score", "Type")
      }
    } else {
      # score_data is already a dataframe
      score_df <- as.data.frame(score_data, stringsAsFactors = FALSE)
      if (!"Type" %in% colnames(score_df)) {
        score_df$Type <- type_labels[i]
      }
      # Ensure first column is named "Score"
      if (colnames(score_df)[1] != "Score") {
        colnames(score_df)[1] <- "Score"
      }
    }


    # Add Genotype_ID column from rownames to preserve it after rbind
    score_df$Genotype_ID <- rownames(score_df)

    # Store in the result list
    result_list[[i]] <- score_df
  }

  # Combine all results - use rbind with make.row.names=FALSE to handle duplicates
  final_result <- do.call(rbind, c(result_list, make.row.names = FALSE))

  return(final_result)
}

#############################################################################

combine_traits_with_groups <- function(processed_data_list) {
  combined_results <- list()

  # Loop through each dataset type
  for (dataset_name in names(processed_data_list)) {
    dataset <- processed_data_list[[dataset_name]]

    # Combine all sublists within the dataset into a single vector
    traits <- unlist(dataset, use.names = FALSE)

    env <- seq(1, length(traits))


    # Create a data frame for this dataset
    combined_results[[dataset_name]] <- data.frame(
      Trait = traits,
      Env = env
    )
  }

  # Combine all datasets into one large data frame
  final_result <- do.call(rbind, combined_results)

  # Add a dataset identifier to distinguish between datasets
  final_result$Dataset <- rep(names(combined_results),
    times = sapply(combined_results, nrow)
  )

  return(final_result)
}


##############################################################################

# Now call preprocess_data with the indices parameter.
processed_data_list <- list(
  gaussian = preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  sinusoidal = preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  wave = preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices),
  linear = preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
)

# Print processing information
total_entries <- sum(sapply(processed_data_list, length))
entries_per_form <- sapply(processed_data_list, length)
message("\n=== Data Processing Summary ===")
message("KEEP_REPLICATES: ", KEEP_REPLICATES)
message("Total genotype entries to process: ", total_entries)
message("  Linear: ", entries_per_form["linear"], " entries")
message("  Gaussian: ", entries_per_form["gaussian"], " entries")
message("  Sinusoidal: ", entries_per_form["sinusoidal"], " entries")
message("  Wave: ", entries_per_form["wave"], " entries")
if (KEEP_REPLICATES) {
  message("(80 genotypes × 3 replicates per form = 240 entries/form)")
} else {
  message("(80 genotypes averaged across replicates = 80 entries/form)")
}
message("Calculating plasticity scores...\n")

linear <- preprocess_data(linear_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
gaussian <- preprocess_data(gaussian_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
sinusoidal <- preprocess_data(sinusoidal_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
wave <- preprocess_data(wave_data, trait_col = 4, genotype_col = 1, replicate_col = 2, indices = indices)
type_labels <- c("gaussian", "sinusoidal", "wave", "linear")


gaussian_long <- combine_traits_with_groups(gaussian)
wave_long <- combine_traits_with_groups(wave)
sinus_long <- combine_traits_with_groups(sinusoidal)
linear_long <- combine_traits_with_groups(linear)

# test <- data.frame(genotype = rep(1, 40),
#                   replicate = c(rep(0, 10), rep(1, 10), rep(2, 10), rep(3, 10)),
#                   Trait = c(rep(1, 20), rep(3, 20)))
# test_processed <- preprocess_data(test, trait_col = 3, genotype_col = 1, replicate_col = 2, indices = indices)
#
# CV_t_test <- call_function(test_processed, calculate_CVt)
# RN_test <- call_function(test_processed, calculate_reaction_norm_slope)
#################################################################################################################################### 1.R
####################################################################################################################################

CV_t <- post_process(processed_data_list, calculate_CVt, type_labels)

RN <- post_process(processed_data_list, calculate_reaction_norm_slope, type_labels)

RNN <- post_process(processed_data_list, calculate_reaction_norm_non_linear, type_labels, 3)

D_slope <- post_process(processed_data_list, calculate_D_slope, type_labels)

RC <- post_process(processed_data_list, calculate_RC, type_labels)

# comparative score -- not one score for a single genotype but for a population of different genotypes. Therefore not suitable for direct comparison
CVm_gaussian <- calculate_CVm(trait_values = gaussian)
CVm_wave <- calculate_CVm(trait_values = wave)
CVm_sinus <- calculate_CVm(trait_values = sinusoidal)
CVm_linear <- calculate_CVm(trait_values = linear)

# comparative score -- not one score for a single genotype but for a population of different genotypes. Therefore not suitable for direct comparison
CVmd_gaussian <- calculate_CVmd(trait_values = gaussian)
CVmd_wave <- calculate_CVmd(trait_values = wave)
CVmd_sinus <- calculate_CVmd(trait_values = sinusoidal)
CVmd_linear <- calculate_CVmd(trait_values = linear)

Covariate <- rep(1, length(indices)) # static covariate to minimize influence on the realisation of the score
# sequential environment


gPi <- post_process(processed_data_list, calculate_grand_plasticity, type_labels, env, Covariate)

PPF <- post_process(processed_data_list, calculate_PPF, type_labels, env, Covariate)

PPi <- post_process(processed_data_list, calculate_Phenotypic_Plasticity_Index, type_labels)

PImd <- post_process(processed_data_list, calculate_PImd, type_labels, env)

PILSM <- post_process(processed_data_list, calculate_PILSM, type_labels, env, Covariate)

RTR <- post_process(processed_data_list, calculate_RTR, type_labels, env)

PIR <- post_process(processed_data_list, calculate_PIR, type_labels, env)


############### 2.R
###############

RDPI <- post_process(processed_data_list, calculate_rdpi, type_labels)

ESPI <- post_process(processed_data_list, calculate_ESPI, type_labels)

ESPIID <- post_process(processed_data_list, calculate_espiid, type_labels)


################# 3.R
#################


PSI <- post_process(processed_data_list, calculate_PSI, type_labels)

RPI <- post_process(processed_data_list, calculate_RPI, type_labels)

PQ <- post_process(processed_data_list, calculate_PQ, type_labels)

PR <- post_process(processed_data_list, calculate_PR, type_labels) # not really comparable as it required individual level data

NRW <- post_process(processed_data_list, calculate_NRW, type_labels)

ESP <- post_process(processed_data_list, calculate_ESP, type_labels)

PD <- post_process(processed_data_list, calculate_general_PD, type_labels) # not comparable as it required stress env and control env

FPI <- post_process(processed_data_list, calculate_FPI, type_labels) # not comparable as it required stress env and control env

################## 4.R
##################


add_arg <- rep(1, length(linear[1][[1]]))

# DPI=post_process(processed_data_list,calculate_DPI,type_labels,add_arg) # not comparable as this requires time resolved data

CEV <- post_process(processed_data_list, calculate_CEV, type_labels)

# PRI=post_process(processed_data_list,calculate_PRI,type_labels) not comparable as extreme env needs to be specified

PFI <- post_process(processed_data_list, calculate_PFI, type_labels)

APC <- post_process(processed_data_list, calculate_APC, type_labels)

SI <- post_process(processed_data_list, calculate_SI, type_labels)

RSI <- post_process(processed_data_list, calculate_RSI, type_labels)

EVS <- post_process(processed_data_list, calculate_EVS, type_labels)

MVPi <- post_process(processed_data_list, calculate_MVPi, type_labels)

################### 5.R
###################

# Plasticity=post_process(processed_data_list,calculate_plasticity,type_labels)
#
# env_cov=post_process(processed_data_list,cross_env_cov,type_labels)
#
scores_list <- list(
  CV_t = CV_t,
  RN = RN,
  RNN = RNN,
  D_slope = D_slope,
  RC = RC,
  gPi = gPi,
  PPF = PPF,
  PPi = PPi,
  PImd = PImd,
  PILSM = PILSM,
  RTR = RTR,
  PIR = PIR,
  RDPI = RDPI,
  ESPI = ESPI,
  ESPIID = ESPIID,
  PSI = PSI,
  RPI = RPI,
  PQ = PQ,
  PR = PR,
  NRW = NRW,
  ESP = ESP,
  CEV = CEV,
  PFI = PFI,
  APC = APC,
  SI = SI,
  RSI = RSI,
  EVS = EVS
)

keys <- names(processed_data_list)
fw_matrices <- setNames(vector("list", length(keys)), keys)

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  L <- max(lengths(gl))
  env_ids <- seq_len(L)
  M <- matrix(NA_real_,
    nrow = length(G), ncol = L,
    dimnames = list(G, paste0("E", env_ids))
  )
  for (g in G) {
    v <- as.numeric(gl[[g]])
    M[g, seq_along(v)] <- v
  }
  fw_matrices[[k]] <- as.data.frame(M)
}

linear_matrix <- fw_matrices$linear
gaussian_matrix <- fw_matrices$gaussian
sinusoidal_matrix <- fw_matrices$sinusoidal
wave_matrix <- fw_matrices$wave

# --- Finlay-Wilkinson (population-level: requires full genotype x environment matrix) ---
message("Calculating Finlay-Wilkinson scores...")
fw_avg_matrices <- lapply(fw_matrices, function(mat) {
  mat <- as.matrix(mat)
  if (exists("KEEP_REPLICATES") && KEEP_REPLICATES) {
    base_names <- sub("_Rep_.*", "", rownames(mat))
    unique_genos <- unique(base_names)
    avg_mat <- do.call(rbind, lapply(unique_genos, function(g) {
      colMeans(mat[base_names == g, , drop = FALSE], na.rm = TRUE)
    }))
    rownames(avg_mat) <- unique_genos
    avg_mat
  } else {
    mat
  }
})

# Per-form FW: each form uses its own population as the reference (used for GWAS and per-form figures)
# FW always produces exactly one beta per genotype (computed on replicate-averaged data).
# We never duplicate that value across replicates — doing so would inflate sample size without
# adding information and is statistically wrong for GWAS.
FW <- do.call(rbind, lapply(names(fw_avg_matrices), function(form) {
  fw_res <- calculate_finlay_wilkinson(fw_avg_matrices[[form]])
  data.frame(
    Score = fw_res$beta,
    Type = form,
    Genotype_ID = fw_res$genotype,
    stringsAsFactors = FALSE
  )
}))

# Cross-form FW: all forms combined into one matrix, single shared environment index
# (used for figures that compare across genotype forms)
# Prefix form name to rownames so all 200 IDs are unique (each form uses genotype IDs 1-50)
combined_fw_mat <- do.call(rbind, lapply(names(fw_avg_matrices), function(form) {
  mat <- fw_avg_matrices[[form]]
  rownames(mat) <- paste0(form, "_", rownames(mat))
  mat
}))
combined_fw_res <- calculate_finlay_wilkinson(combined_fw_mat)
form_labels <- rep(names(fw_avg_matrices), times = sapply(fw_avg_matrices, nrow))
FW_combined <- data.frame(
  Score = combined_fw_res$beta,
  Type = form_labels,
  Genotype_ID = combined_fw_res$genotype,
  stringsAsFactors = FALSE
)

fw_long_all <- data.frame(
  Dataset = character(), genotype = character(),
  environment = integer(), y = double(),
  stringsAsFactors = FALSE
)

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  for (g in G) {
    v <- as.numeric(gl[[g]])
    n <- length(v)
    df <- data.frame(Dataset = k, genotype = g, environment = seq_len(n), y = v)
    fw_long_all <- rbind(fw_long_all, df)
  }
}


keys <- names(processed_data_list)
row_list <- list()
fw_long_all <- data.frame(Dataset = character(), genotype = character(), environment = integer(), y = double(), stringsAsFactors = FALSE)
Lmax <- max(unlist(lapply(processed_data_list, function(gl) max(lengths(gl)))))

for (k in keys) {
  gl <- processed_data_list[[k]]
  G <- names(gl)
  for (g in G) {
    uid <- paste(k, g, sep = "_")
    v <- as.numeric(gl[[g]])
    n <- length(v)
    row_list[[uid]] <- c(v, rep(NA_real_, Lmax - n))
    fw_long_all <- rbind(fw_long_all, data.frame(Dataset = k, genotype = uid, environment = seq_len(n), y = v))
  }
}


scores_list$FW <- FW

# -----------------------------------------------------------------------------
# SAVE ALL DATA TO DISK (Timestamped)
# -----------------------------------------------------------------------------
output_dir <- paste(OUTPUT_BASE,"/synthetic_data/scores_output",sep="")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
run_dir <- file.path(output_dir, paste0("run_", timestamp))
dir.create(run_dir)

message("Saving all simulation data to: ", run_dir)

# 1. Scores
saveRDS(scores_list, file.path(run_dir, "scores_list.rds"))

# 2. Reaction Norms (Phenotypes)
saveRDS(combined_df, file.path(run_dir, "reaction_norms.rds"))
write.csv(combined_df, file.path(run_dir, "reaction_norms.csv"), row.names = FALSE)

# 3. Genetics (Copy from temp location or save directly if variables exist)
if (exists("snp_matrix")) {
  saveRDS(snp_matrix, file.path(run_dir, "snp_matrix.rds"))
  write.csv(snp_matrix, file.path(run_dir, "snp_matrix.csv"))
}
if (exists("params_df")) {
  write.csv(params_df, file.path(run_dir, "genotypic_parameters.csv"), row.names = FALSE)
}
if (exists("truth_df")) {
  write.csv(truth_df, file.path(run_dir, "causal_snps_truth.csv"), row.names = FALSE)
}

# 4. Kinship Matrix (Calculate if not exists)
if (!exists("K")) {
  if (exists("snp_matrix")) {
    message("Calculating Kinship Matrix...")
    # Simple centered IBS or use sommer if available
    if (requireNamespace("sommer", quietly = TRUE)) {
      # sommer expects -1, 0, 1 usually, or 0, 1, 2. Our data is 0, 1, 2.
      # A.mat expects -1, 0, 1. Let's subtract 1.
      M_adj <- snp_matrix - 1
      K <- sommer::A.mat(M_adj)
    } else {
      # Fallback: simple correlation or identity
      K <- cor(t(snp_matrix))
    }
  }
}

if (exists("K")) {
  saveRDS(K, file.path(run_dir, "kinship_matrix.rds"))
}

message("Done. All data saved in: ", run_dir)

# -----------------------------------------------------------------------------
# SAVE PER-PPI SCORE CSVs  (one file per score, genotype → trait)
# -----------------------------------------------------------------------------
ppi_dir <- file.path(run_dir, "ppi_scores")
dir.create(ppi_dir, showWarnings = FALSE)

for (score_name in names(scores_list)) {
  df <- scores_list[[score_name]]
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!"Genotype_ID" %in% colnames(df)) df$Genotype_ID <- rownames(df)
  out <- df[, c("Genotype_ID", "Score", "Type"), drop = FALSE]
  colnames(out)[colnames(out) == "Score"] <- score_name
  write.csv(out, file.path(ppi_dir, paste0(score_name, ".csv")), row.names = FALSE)
}
message("Per-PPI CSVs saved to: ", ppi_dir)

# NOTE on FW: FW produces exactly one value per genotype (averaged over replicates).
# Its CSV will therefore have fewer rows than other scores when KEEP_REPLICATES=TRUE.
# This is correct — never use #replicates copies of the same FW beta in GWAS.

# -----------------------------------------------------------------------------
# SAVE PLINK-COMPATIBLE FILES  (.map + .ped)
# -----------------------------------------------------------------------------
if (exists("snp_matrix")) {
  plink_dir <- file.path(OUTPUT_BASE, "synthetic_data", "genetics")

  # .map  — Chr1  SNP_i  0  position(bp)
  # Positions on a fixed raster of 345 bp so tools that care about spacing work sensibly.
  snp_ids  <- colnames(snp_matrix)
  snp_pos  <- seq_along(snp_ids) * 345L
  map_df <- data.frame(
    chr = "Chr1",
    id  = snp_ids,
    cm  = 0,
    pos = snp_pos,
    stringsAsFactors = FALSE
  )
  write.table(map_df,
              file.path(plink_dir, "snp_data.map"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

  # .ped  — FID IID PaternalID MaternalID Sex Phenotype [allele pairs per SNP]
  # Allele encoding (raw snp_matrix values 0/1/2, i.e. before rrBLUP's -1 shift):
  #   0  (hom ref,  -1 in rrBLUP)  → A A
  #   1  (het,       0 in rrBLUP)  → A C
  #   2  (hom alt,  +1 in rrBLUP)  → C C
  allele_map <- c("0" = "A A", "1" = "A C", "2" = "C C")

  ped_lines <- character(nrow(snp_matrix))
  geno_names <- rownames(snp_matrix)
  for (i in seq_len(nrow(snp_matrix))) {
    alleles <- allele_map[as.character(snp_matrix[i, ])]
    alleles[is.na(alleles)] <- "0 0"          # encode any unexpected value as missing
    ped_lines[i] <- paste(
      geno_names[i], geno_names[i], "0", "0", "0", "-9",
      paste(alleles, collapse = "\t"),
      sep = "\t"
    )
  }
  writeLines(ped_lines, file.path(plink_dir, "snp_data.ped"))

  message("PLINK files saved to: ", plink_dir)
  message("  snp_data.map  (", nrow(map_df), " SNPs)")
  message("  snp_data.ped  (", nrow(snp_matrix), " genotypes)")
}
