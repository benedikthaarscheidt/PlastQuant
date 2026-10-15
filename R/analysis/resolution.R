# =============================================================================
# resolution.R — resolution experiment (dataprep + regression combined)
# =============================================================================
# WHAT IT DOES: Runs, in one pass, what regression_sumstats_dataprep.R and
#   regression_summarystats.R do together, for the resolution experiment: Kendall
#   correlations (permutation test) between every score and every summary statistic,
#   pooled over the forms and per form, for every sub-range and sampling interval.
# REQUIRES:     POPULATION_DIR: a simulated population folder; its reaction norms are
#               sampled at every range x interval of a grid of global_initial_length points.
# PRODUCES:     correlations_<range>_interval_<i>_indices_<...>.csv in
#               POPULATION_DIR/analyses/correlation_summary_stats/ (read by plotting_resolution.R).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "analysis", "resolution.R"))
# -----------------------------------------------------------------------------
# PARAMETERS: inherits 03's parameters via the auto-source; resolution knobs are set in
#   the body of the script (edit them there).
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

# Purpose:   This script does what regression_sumstats_dataprep.R and regression_summarystats.R do in combination).
#            For each specified sub‐range of the trait curve and each sampling
#            interval:
#            1. Sample the population's reaction norms at the new resolution and
#               sub‐range (sample_population(); same genotypes and plants),
#            2. Compute the plasticity scores (compute_ppis()),
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
#   • POPULATION_DIR: the simulated population folder
#   • R/pipeline: sampling, PPIs, summary_statistics(), calc_similarity_perm()
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

# Define summary statistic names (produced by summary_statistics())
stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")

# The population whose reaction norms are sampled at every design
if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a simulated population folder) first")
source(here::here("R", "pipeline", "load_pipeline.R"))
correlation_dir <- file.path(POPULATION_DIR, "analyses", "correlation_summary_stats")
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

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
    indices <- grid_design(global_start_index, global_end_index, global_sampling_interval)
    
    #print(indices)
    
    # The population's reaction norms sampled at this design (same genotypes and plants)
    measured <- sample_population(POPULATION_DIR, t = design_points(indices, global_initial_length))
    genotype_forms <- c("linear", "gaussian", "sinusoidal", "wave")

    # The 27 plasticity scores of this analysis (all PPIs except the population-level FW)
    scores_list <- compute_ppis(measured$phenotypes, measured$environments, exclude = "FW")

    # Replicate-mean reaction norm matrices of every genotype form (rows in the order of the scores)
    rn_matrix_list <- lapply(setNames(genotype_forms, genotype_forms), function(form) {
      replicate_means(measured$phenotypes[[form]], measured$environments)
    })

    # Helper: extract score vector for a given genotype form.
    prepare_score_vector <- function(score_df, genotype_form) {
      subset(score_df, score_df[, "Type"] == genotype_form)[, "Score"]
    }

    # For global analysis, combine summary stats from all genotype forms.
    combined_stats <- data.frame()
    for (form in genotype_forms) {
      form_stats <- summary_statistics(rn_matrix_list[[form]], measured$environments)
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
        stats_df <- summary_statistics(rn_matrix_list[[form]], measured$environments)
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
    csv_filename <- file.path(correlation_dir, paste0("correlations_", range_name, "_interval_", interval, "_indices_", paste(indices, collapse = "_"), ".csv"))
    write.csv(all_results, file = csv_filename, row.names = FALSE)
    message("CSV saved as ", csv_filename)
    csv_files <- c(csv_files, csv_filename)
    
    cat("  Completed analysis for sampling interval =", interval, "with range", range_name, "\n\n")
  }
}

cat("All analyses completed. CSV files saved:\n")
print(csv_files)
