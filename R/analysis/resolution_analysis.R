# =============================================================================
# resolution_analysis.R — sampling sensitivity: correlation and its regression
# =============================================================================
# WHAT IT DOES: From the tables of regression_sumstats_dataprep.R (the same reaction
#   norms sampled at 2, 3, 4, 5, 6, 11, 26 and 50 points) computes
#     - the Pearson correlation of every PPI with every summary statistic across all
#       genotypes, per number of sampling points,
#     - for every PPI-statistic pair the linear regression of that correlation on the
#       number of sampling points (slope, its standard error, t-test).
#   The figure script R/figures/combined_resolution_figure.R reads these results and
#   fits the mixed model on the tables.
# REQUIRES:     regression_sumstats_dataprep.R run on POPULATION_DIR.
# PRODUCES:     POPULATION_DIR/analyses/resolution/
#                 correlations.csv            N_points, PPI, Stat, r
#                 correlation_regression.csv  PPI, Stat, intercept, slope, se_slope, t_stat, p_value, nobs
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R", "analysis", "resolution_analysis.R"))
# -----------------------------------------------------------------------------
# PARAMETERS
#   POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder               [COMMON]
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")

tables <- read_regression_tables(file.path(POPULATION_DIR, "analyses", "regression_summary_stats"), range = "full")
out_dir <- population_subdir(POPULATION_DIR, "analyses", "resolution")
cors <- correlations_by_resolution(tables)
write.csv(cors, file.path(out_dir, "correlations.csv"), row.names = FALSE)
write.csv(correlation_regression(cors), file.path(out_dir, "correlation_regression.csv"), row.names = FALSE)
message("Resolution analysis written to ", out_dir)
