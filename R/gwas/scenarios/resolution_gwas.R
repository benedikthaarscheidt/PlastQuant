# =============================================================================
# resolution_gwas.R — influence of the resolution on the GWAS (command-line steps)
# =============================================================================
# WHAT IT DOES: Steps of the resolution GWAS, called by hpc/resolution_gwas.slurm (or by
#   hand). The resolutions themselves come from regression_sumstats_dataprep.R, so the
#   GWAS runs on the same resolution folders as the regression analyses.
#     simulate  simulates the population into <population folder> (data/, variance/);
#               an existing population is kept
#     gwas      GWAS of one resolution folder n_<n> with its power
#               (n_<n>/<handling>/gwas/); computes the PPIs of the replicate handling
#               first if the folder has none
#   Every GWAS only reads data/ and writes its own n_<n>/ folder, so the resolutions can
#   run at the same time. combine_resolution_gwas.R collects the power afterwards.
# HOW TO RUN:   from the repo root
#   Rscript --vanilla R/gwas/scenarios/resolution_gwas.R simulate <population folder> <h2> <Mc> <seed>
#   Rscript --vanilla R/gwas/scenarios/resolution_gwas.R gwas <resolution folder> [mean_before|replicates]
# -----------------------------------------------------------------------------
# ARGUMENTS
#   h2        numeric in [0, 1]   share of the latent genetic variance from the causal SNPs
#   Mc        integer             causal SNPs per parameter
#   seed      integer             seed of the population
#   handling  "mean_before" (default) or "replicates"   replicate handling of the PPIs
# All other settings are the defaults of sim_config() in R/sim/simulation_config.R.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

args <- commandArgs(trailingOnly = TRUE)
usage <- paste("usage: resolution_gwas.R simulate <population folder> <h2> <Mc> <seed>\n",
               "      resolution_gwas.R gwas <resolution folder> [mean_before|replicates]")
if (length(args) < 2 || !args[1] %in% c("simulate", "gwas")) stop(usage)
source(here::here("R", "pipeline", "load_pipeline.R"))
t0 <- Sys.time()

if (args[1] == "simulate") {
  if (length(args) < 5) stop(usage)
  population_dir <- args[2]
  h2 <- as.numeric(args[3]); Mc <- as.integer(args[4]); seed <- as.integer(args[5])
  if (anyNA(c(h2, Mc, seed))) stop(usage)
  if (file.exists(file.path(population_dir, "data", "norms.rds"))) {
    message("Population exists, kept: ", population_dir)
  } else {
    simulate_population(sim_config(h2 = h2, Mc = Mc, seed = seed), population_dir)   # -> data/, variance/
    message("Population simulated in ", format(round(Sys.time() - t0, 1)), ": ", population_dir)
  }
} else {
  resolution_dir <- args[2]
  replicate_handling <- if (length(args) >= 3) args[3] else "mean_before"
  if (!file.exists(file.path(ppi_set_path(resolution_dir, replicate_handling), "ppis", "ppi_scores.rds"))) {
    compute_population_ppis(resolution_dir, replicate_handling)               # -> <handling>/ppis/
  }
  run_population_gwas(resolution_dir, replicate_handling)                     # -> <handling>/gwas/, plots/
  message("GWAS of ", basename(resolution_dir), " done in ", format(round(Sys.time() - t0, 1)))
}
