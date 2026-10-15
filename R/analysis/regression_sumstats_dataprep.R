# =============================================================================
# regression_sumstats_dataprep.R — PPIs of one population at all resolutions
# =============================================================================
# WHAT IT DOES: For every resolution n in RESOLUTIONS measures the true reaction norms of
#   one simulated population at n equidistant points (same genotypes and plants; only
#   the sampling changes) and computes all plasticity scores and the nine summary
#   statistics. Every resolution is written as its own folder, exactly as when done one
#   by one (sample_resolution(), compute_population_ppis()). A resolution that already
#   has its PPIs is kept. These folders are the input of resolution_analysis.R,
#   combined_resolution_figure.R, correlation_resolution.R and
#   correlation_resolution_experiment.R.
# REQUIRES:     POPULATION_DIR: a simulated population folder (simulate_population()).
# PRODUCES:     POPULATION_DIR/n_<n>/ for every n: data/ (measurements), variance/,
#               ppis/ (scores), analyses/summary_statistics.csv
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R","analysis","regression_sumstats_dataprep.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — POPULATION_DIR is set before sourcing; the resolutions are edited here
#   POPULATION_DIR  path, REQUIRED   simulated population folder                         [COMMON]
#   RESOLUTIONS     integer, default c(2, 3, 4, 5, 6, 11, 26, 50)  sampling points        [COMMON]
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a simulated population folder) first")
RESOLUTIONS <- c(2, 3, 4, 5, 6, 11, 26, 50)

res_dirs <- compute_resolutions(POPULATION_DIR, RESOLUTIONS)
message("PPIs and summary statistics of ", length(res_dirs), " resolutions in ", POPULATION_DIR)
