# =============================================================================
# regression_sumstats_dataprep.R — build the summary-stats vs scores tables
# =============================================================================
# WHAT IT DOES: For every sampling design (range of the gradient x sampling interval on
#   a grid of global_initial_length points) samples the reaction norms of one simulated
#   population again (same genotypes and plants; only the sampling points change),
#   computes all plasticity scores and the nine summary statistics per genotype, and
#   writes one table per design. These tables are the input of
#   resolution_analysis.R, regression_summarystats.R, kendall_rank_corr.R,
#   kendall_tau_corr_clustering.R and the resolution figures.
# REQUIRES:     POPULATION_DIR: a simulated population folder (simulate_population()).
# PRODUCES:     POPULATION_DIR/analyses/regression_summary_stats/
#                 regression_data_<range>_interval_<interval>_indices_<indices>.csv
#               (first column the genotype, then one column per score, then the
#               summary statistics min ... mean_upper).
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R","analysis","regression_sumstats_dataprep.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — POPULATION_DIR is set before sourcing; the design settings are edited here
#   POPULATION_DIR         path, REQUIRED   simulated population folder                 [COMMON]
#   sampling_intervals     integer, default c(1, 2, 5, 10, 15, 20, 25, 49)
#                          (2, 3, 4, 5, 6, 11, 26 and 50 sampling points)               [COMMON]
#   global_initial_length  integer, default 50       grid the indices refer to
#   range_list             list, default list(full = c(1, 50))  ranges of the grid
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a simulated population folder) first")
sampling_intervals <- c(1, 2, 5, 10, 15, 20, 25, 49)
global_initial_length <- 50
range_list <- list(full = c(1, global_initial_length))

output_folder <- file.path(POPULATION_DIR, "analyses", "regression_summary_stats")
files <- write_regression_data(POPULATION_DIR, range_list, sampling_intervals, output_folder,
                               grid_length = global_initial_length)
message("Regression data saved to ", output_folder, " (", length(files), " tables)")
