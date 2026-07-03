# =============================================================================
# regression_sumstats_dataprep.R — build the summary-stats vs scores table
# =============================================================================
# WHAT IT DOES: Assembles the table of reaction-norm summary statistics against
#   plasticity scores that regression_summarystats.R and several figure scripts consume.
#   Auto-sources 03_plasticity_scores.R and analysis/summary_stats_vs_scores.R.
# REQUIRES:     Runs 03 (so it needs HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE set
#               by a scenario or your session, as 03 does).
# PRODUCES:     Regression-input CSV(s) under output/regression_summary_stats/.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","regression_sumstats_dataprep.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (supplied by 03 / a scenario / your session; see 03_plasticity_scores.R)
#   NUM_GENOTYPES, HERITABILITY_5TH, CAUSAL_SNP_NUM, USE_GENETICS, OUTPUT_BASE      [COMMON]
# =============================================================================


# Script:    regression_sumstats_dataprep.R
# Purpose:   For each specified sampling interval (resolution) and range of the
#            reaction‐norm curve, this script:
#              1. Defines the sampling indices over the full or sub‐range;
#              2. Sources head3.R and summary_stats_vs_scores.R to recalc reaction‑
#                 norm matrices and plasticity scores at that resolution;
#              3. Computes nine summary statistics (min, max, mean, … mean_upper)
#                 for each genotype form (linear, gaussian, sinusoidal, wave);
#              4. Builds the predictor matrix X (summary stats) and the response
#                 matrix y (one column per plasticity score);
#              5. Combines X and y into one data frame and writes it to CSV.
#
# Inputs & Globals:
#   • sampling_intervals: vector of step‐sizes (1,2,5,…,49)
#   • global_initial_length: length of full trait vector (50)
#   • range_list: named list of (start, end) index pairs (here just “full”)
#   • head3.R: re‐samples data & computes score_list, rn_matrix_list, genotype_forms
#   • summary_stats_vs_scores.R: provides compute_summary_stats() and prepare_score_vector()
#
# Dependencies:
#   factoextra, dbscan, cluster, mclust, pheatmap, corrplot, vegan, energy,
#   gridExtra, grid, gtable
#
# Outputs:
#   • One CSV per (range × interval) in output_folder, named
#       regression_data_<range>_interval_<interval>_indices_<indices>.csv
#     containing columns = all score values + all summary‐stat predictors.

# ---- Genetics configuration -----------------------------------------------
# Set USE_GENETICS = FALSE to skip SNP simulation entirely.
# Reaction norms are then generated with purely random parameters (no causal
# variants, no heritability). Set to TRUE to run the full genetic pipeline
# (requires HERITABILITY_5TH, CAUSAL_SNP_NUM, and OUTPUT_BASE).
USE_GENETICS  <- FALSE
OUTPUT_BASE   <- here::here("output", "no_genetics_output")
NUM_GENOTYPES <- 200   # must be divisible by 4; keep small for fast runs
# ---------------------------------------------------------------------------

ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

output_folder <- here::here("output", "regression_summary_stats")
ensure_dir_exists(output_folder)

sampling_intervals <- c(1, 2, 5, 10, 15, 20,25, 49)
global_initial_length <- 50


range_list <- list(
  full = c(1, global_initial_length)
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





stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")

for (range_name in names(range_list)) {
  global_start_index <- range_list[[range_name]][1]
  global_end_index   <- range_list[[range_name]][2]
  
  cat("Processing range:", range_name, "(", global_start_index, "to", global_end_index, ")\n")
  
  # Loop over each sampling interval (each resolution)
  for (interval in sampling_intervals) {
    cat("  Processing resolution (sampling interval):", interval, "\n")
    
    # Determine the indices used to subset the data (ensure the end is included)
    indices <- unique(c(seq(global_start_index, global_end_index, by = interval),
                        global_end_index))
    env = seq(1, 10,length.out=length(indices))
    traits <- length(indices)
    
    
    source(here::here("R", "03_plasticity_scores.R"))
    data_loaded <- TRUE
    source(here::here("R", "analysis", "summary_stats_vs_scores.R"))
    data_loaded <- FALSE  # reset so head3.R re-runs on the next interval iteration
    
    ###################################
    ## Build Predictor Matrix (X)
    ###################################
    # For each genotype, compute its summary statistics on the subset of data defined by indices.
    X <- data.frame()
    for (gen in genotype_forms) {
      rn_matrix <- rn_matrix_list[[gen]]
      print(gen)
      
      stats <- compute_summary_stats(rn_matrix)[, stat_names, drop = FALSE]
    
      stats$genotype <- gen
      
      X <- rbind(X, stats)
    }
    
    
    ###################################
    ## Build the Response Matrix (y)
    ###################################
   
    num_gen <- nrow(X)
    num_scores <- length(scores_list)
    y_matrix <- matrix(NA, nrow = num_gen, ncol = num_scores)
    print(dim(y_matrix))
    colnames(y_matrix) <- names(scores_list)
    rownames(y_matrix) <- seq(1,nrow(X))
    
    
    
    for(s in seq_along(scores_list)){
      score_name <- names(scores_list)[s]
      score_df <- scores_list[[score_name]]
      score_value=c()
      for (i in seq_along(genotype_forms)) {
        gen <- genotype_forms[i]
        score_value <- c(score_value,as.numeric(as.character(prepare_score_vector(score_df, gen))))
      }
      
      y_matrix[,s ] <- score_value
      
    }
    
    ###################################
    ## Save the Regression Data
    ###################################
    
    regression_data <- list(
      X_full = X[, stat_names, drop = FALSE],
      y = y_matrix
    )
    combined_df <- cbind(regression_data$y,regression_data$X_full)
    
    save_filename <- file.path(output_folder, paste0("regression_data_", range_name, "_interval_", interval, "_indices_", paste(indices, collapse = "_" ),".csv"))
    write.csv(combined_df, file = save_filename, row.names = TRUE,col.names = TRUE)
    message("Regression data saved to ", save_filename)
    
  } # end for each sampling interval
} # end for each range

cat("All regression data have been computed and saved.\n")
