# =============================================================================
# resolution_analysis.R — sampling sensitivity of the PPIs (computation)
# =============================================================================
# WHAT IT DOES: Samples the reaction norms of one simulated population again at every
#   resolution (same genotypes and plants; only the sampling points change), computes
#   the 28 PPIs and the nine summary statistics per genotype, and from these
#     - the Pearson correlation of every PPI with every summary statistic across all
#       genotypes, per resolution,
#     - the linear regression of that correlation on the number of sampling points.
#   The figure script (R/figures/combined_resolution_figure.R) reads these results and
#   fits the mixed model on the tables.
# REQUIRES:     POPULATION_DIR: a simulated population folder (simulate_population()).
# PRODUCES:     POPULATION_DIR/analyses/resolution/
#                 tables/n_<points>.csv     PPIs and summary statistics per genotype
#                 correlations.csv          N_points, PPI, Stat, r
#                 correlation_regression.csv  PPI, Stat, slope, se_slope, t_stat, p_value, ...
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R", "analysis", "resolution_analysis.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   POPULATION_DIR  path, REQUIRED                              simulated population      [COMMON]
#   RESOLUTIONS     integer, default c(2, 3, 4, 5, 6, 11, 26, 50)  sampling points per design [COMMON]
#   EXCLUDE_PPIS    character, default NULL                      PPIs left out
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a simulated population folder) first")
if (!exists("RESOLUTIONS")) RESOLUTIONS <- c(2, 3, 4, 5, 6, 11, 26, 50)
if (!exists("EXCLUDE_PPIS")) EXCLUDE_PPIS <- NULL

registry <- default_ppi_registry()
registry <- registry[setdiff(names(registry), EXCLUDE_PPIS)]

out_dir <- file.path(POPULATION_DIR, "analyses", "resolution")
unlink(out_dir, recursive = TRUE)
table_dir <- population_subdir(POPULATION_DIR, "analyses", "resolution", "tables")

tables <- resolution_tables(POPULATION_DIR, RESOLUTIONS, registry)
for (tab in tables) {
  write.csv(tab, file.path(table_dir, paste0("n_", tab$N_points[1], ".csv")), row.names = FALSE)
}
cors <- correlations_by_resolution(tables)
write.csv(cors, file.path(out_dir, "correlations.csv"), row.names = FALSE)
write.csv(correlation_regression(cors), file.path(out_dir, "correlation_regression.csv"), row.names = FALSE)
message("Resolution analysis written to ", out_dir)
