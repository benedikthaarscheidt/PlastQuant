# =============================================================================
# combine_scenarios.R — collect the heritability scenarios into one table set
# =============================================================================
# WHAT IT DOES: Takes the population folder of every scenario_<h2>_<Mc>, keeps the
#   GWAS hits below P_THRESHOLD and the causal SNPs (both tagged with H2 and Mc), and
#   counts the found causal SNPs per scenario, shape, parameter and PPI. Every shape
#   is evaluated against its own causal SNPs only (see detection_counts()).
# REQUIRES:     Scenario runs written by R/scenarios/scenario_*.R.
# PRODUCES:     Under OUT_DIR: gwas_hits.csv, causal_snps_truth.csv, tested.csv,
#               detection_counts.csv, populations.csv (read by 06_analyze_gwas_results.R).
# HOW TO RUN:   Rscript --vanilla R/scenarios/combine_scenarios.R   (from the repo root)
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   SCENARIO_ROOT  path, default output/populations/       folder with the scenario_* populations
#   SCENARIOS      character, default all scenario_<h2>_<Mc> folders in SCENARIO_ROOT
#   P_THRESHOLD    numeric, default 1e-3                   P value below which a SNP is a hit [COMMON]
#   N_POINTS       integer, default 6                      resolution n_<N_POINTS> whose GWAS is used
#   OUT_DIR        path, default output/combined_scenarios where the tables are written
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("SCENARIO_ROOT")) SCENARIO_ROOT <- here::here("output", "populations")
if (!exists("SCENARIOS")) {
  SCENARIOS <- grep("^scenario_[0-9]+_[0-9]+$", list.files(SCENARIO_ROOT), value = TRUE)
}
if (!exists("P_THRESHOLD")) P_THRESHOLD <- 1e-3
if (!exists("N_POINTS")) N_POINTS <- 6
if (!exists("OUT_DIR")) OUT_DIR <- here::here("output", "combined_scenarios")
if (!length(SCENARIOS)) stop("no scenario_<h2>_<Mc> folders in ", SCENARIO_ROOT)

message("Combining: ", paste(SCENARIOS, collapse = ", "))
combined <- combine_scenario_runs(file.path(SCENARIO_ROOT, SCENARIOS), n_points = N_POINTS, p_threshold = P_THRESHOLD)
detections <- detection_counts(combined)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(combined$hits, file.path(OUT_DIR, "gwas_hits.csv"), row.names = FALSE)
write.csv(combined$truth, file.path(OUT_DIR, "causal_snps_truth.csv"), row.names = FALSE)
write.csv(combined$tested, file.path(OUT_DIR, "tested.csv"), row.names = FALSE)
write.csv(combined$populations, file.path(OUT_DIR, "populations.csv"), row.names = FALSE)
write.csv(detections, file.path(OUT_DIR, "detection_counts.csv"), row.names = FALSE)
saveRDS(list(p_threshold = P_THRESHOLD, n_points = N_POINTS), file.path(OUT_DIR, "settings.rds"))
message("Tables written to ", OUT_DIR)
