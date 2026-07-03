# Purpose:   This script does what regression_sumstats_dataprep.R and regression_summarystats.R do in combination).
#            For each specified sub‐range of the trait curve and each sampling
#            interval:
#            1. Re‐sample the trait data at the new resolution and sub‐range
#               (via head3.R),
#            2. Re‐compute an array of plasticity scores (via head3.R),
#            3. Compute nine summary statistics on each re‐sampled reaction norm,
#            4. Test Correlation (with 999 permutations) between every score and
#               every summary statistic:
#               • “Global” correlations pooled across all genotype forms,
#               • “Local” correlations within each genotype form (linear,
#                 gaussian, sinusoidal, wave),
#            5. Export a combined table of global+local results to a CSV whose
#               name encodes the range, interval, and actual indices used.
# 
# Inputs:
#   • sampling_intervals: vector of step‐sizes for sampling the curve
#   • global_initial_length: full vector length (default 50)
#   • range_list: named list of (start, end) index pairs
#   • head3.R: re‐sampling & score‐calculation script (uses globals)
#   • summary_stats_vs_scores.R: prepares summary‐stat and correlation functions
#
# Outputs:
#   • One CSV per (range × interval) containing columns:
#       resolution, score, stat, analysis (global/local),
#       genotype_form, correlation (τ), p_value
#   • Printed console log of progress and saved filenames
#
# -----------------------------------------------------------------------
# --- Define sampling intervals and range specifications ---
sampling_intervals <- c ( 49)  # sampling intervals   
global_initial_length <- 50  # original full length of the trait vector

# Define a list of range specifications.
# Each element is a named vector: c(start_index, end_index)
range_list <- list(
  full = c(1, global_initial_length),     # full range: 1 to 50
  partial1 = c(1, global_initial_length/2),    # first half
  partial2 = c(global_initial_length/2, global_initial_length)  # second half
)


library(factoextra)
library(dbscan)
library(cluster)
library(mclust)
library(pheatmap)
library(corrplot)
library(vegan)
library(energy)
library(gridExtra)
library(grid)
library(gtable)

#  Helper function for significance stars
sig_stars <- function(p) {
  if (is.na(p)) {
    return("NA")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Define summary statistic names (produced by compute_summary_stats)
stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")

# Container to track output CSV filenames
csv_files <- c()

# Outer loop: iterate over range specifications
for (range_name in names(range_list)) {
  
  # Set global start and end indices for sampling based on the current range specification
  global_start_index <- range_list[[range_name]][1]
  global_end_index   <- range_list[[range_name]][2]
  
  cat("Running analyses for range:", range_name, "(", global_start_index, "to", global_end_index, ")\n")
  
  # Inner loop: iterate over each sampling interval
  for (interval in sampling_intervals) {
    cat("  Sampling interval =", interval, "\n")
    
    # Set global sampling interval
    global_sampling_interval <- interval
    # Compute the effective indices within the specified range:
    indices <- unique(c(seq(global_start_index, global_end_index, by = global_sampling_interval),
                        global_end_index))
    
    #print(indices)
    
    global_traits <- length(indices)
    traits <- global_traits  
    env=seq_along(indices)
    # Source the processing script that uses the globals.
    # head3.R is modified to use global_sampling_interval, global_start_index, global_end_index, and traits.
    source("~/CRC_1644_Z2_GWAS_simple/R-files/head3.R")
    source("~/CRC_1644_Z2_GWAS_simple/R-files/summary_stats_vs_scores.R")
    # After sourcing, all preprocessed data and score objects (e.g. CV_t, RN, etc.) are recalculated.
    
    # Build the scores list from the computed objects
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
    
    # Combine reaction norm matrices for each genotype form.
    # (Assumes head3.R produces lists: linear, gaussian, sinusoidal, wave)
    linear_matrix   <- do.call(rbind, linear)
    gaussian_matrix <- do.call(rbind, gaussian)
    sinusoidal_matrix <- do.call(rbind, sinusoidal)
    wave_matrix     <- do.call(rbind, wave)
    
    rn_matrix_list <- list(
      linear = linear_matrix,
      gaussian = gaussian_matrix,
      sinusoidal = sinusoidal_matrix,
      wave = wave_matrix
    )
    
    genotype_forms <- c("linear", "gaussian", "sinusoidal", "wave")
    
    # Helper: extract score vector for a given genotype form.
    prepare_score_vector <- function(score_df, genotype_form) {
      subset(score_df, score_df[, "Type"] == genotype_form)[, "Score"]
    }
    
    # For global analysis, combine summary stats from all genotype forms.
    combined_stats <- data.frame()
    for (form in genotype_forms) {
      rn_matrix <- rn_matrix_list[[form]]
      form_stats <- compute_summary_stats(rn_matrix)
      form_stats$Form <- form
      combined_stats <- rbind(combined_stats, form_stats)
    }
    
    # Initialize containers for correlation results.
    # Global analysis: one row per (score, summary stat)
    global_results <- data.frame()
    # Local analysis: one row per (score, summary stat, genotype form)
    local_results <- data.frame()
    
    ## GLOBAL ANALYSIS
    for (score_name in names(scores_list)) {
      score_df <- scores_list[[score_name]]
      overall_scores <- c()
      # Concatenate scores from all genotype forms.
      for (form in genotype_forms) {
        score_subset <- as.numeric(as.character(prepare_score_vector(score_df, form)))
        overall_scores <- c(overall_scores, score_subset)
      }
      for (stat in stat_names) {
        test <- calc_similarity_perm(overall_scores, combined_stats[[stat]], method = "kendall", n_perm = 999)
        global_results <- rbind(global_results, data.frame(
          resolution = interval,
          score = score_name,
          stat = stat,
          analysis = "global",
          genotype_form = "combined",
          correlation = test$estimate,
          p_value = test$p.value
        ))
      }
    }
    
    ## LOCAL ANALYSIS: Compute correlation separately for each genotype form.
    for (score_name in names(scores_list)) {
      score_df <- scores_list[[score_name]]
      for (form in genotype_forms) {
        score_vec <- as.numeric(as.character(prepare_score_vector(score_df, form)))
        rn_matrix <- rn_matrix_list[[form]]
        stats_df <- compute_summary_stats(rn_matrix)
        for (stat in stat_names) {
          test <- calc_similarity_perm(score_vec, stats_df[[stat]], method = "kendall", n_perm = 999)
          local_results <- rbind(local_results, data.frame(
            resolution = interval,
            score = score_name,
            stat = stat,
            analysis = "local",
            genotype_form = form,
            correlation = test$estimate,
            p_value = test$p.value
          ))
        }
      }
    }
    
    # Combine global and local results for this (interval, range) combination.
    all_results <- rbind(global_results, local_results)
    
    # Save the results into a CSV file whose name includes the sampling interval and range label.
    csv_filename <- paste0("correlations_", range_name, "_interval_", interval, "_indices_", paste(indices, collapse = "_"), ".csv")
    write.csv(all_results, file = csv_filename, row.names = FALSE)
    message("CSV saved as ", csv_filename)
    csv_files <- c(csv_files, csv_filename)
    
    cat("  Completed analysis for sampling interval =", interval, "with range", range_name, "\n\n")
  }
}

cat("All analyses completed. CSV files saved:\n")
print(csv_files)
