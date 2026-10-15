# =============================================================================
# resolution_gwas.R — influence of the resolution on the GWAS (command-line driver)
# =============================================================================
# WHAT IT DOES: One population, measured at several resolutions, each with its PPIs and
#   GWAS, so the GWAS power can be compared between resolutions on the SAME true
#   reaction norms. Two steps, called separately (the second one once per resolution,
#   e.g. as a SLURM array, see hpc/):
#     simulate    simulates the population into output/populations/resolution_<h2>_<Mc>_seed<seed>/
#                 (data/, variance/)
#     resolution  measures it at n points (n_<n>/), computes the PPIs
#                 (n_<n>/<handling>/ppis/) and runs the GWAS with its power
#                 (n_<n>/<handling>/gwas/gwas_power.csv)
#   The resolutions only read data/ and write their own n_<n>/ folder, so they can run
#   at the same time. combine_resolution_gwas.R collects the power over the resolutions.
# REQUIRES:     for "resolution": the simulated population (step "simulate")
# HOW TO RUN:   from the repo root
#   Rscript --vanilla R/gwas/scenarios/resolution_gwas.R simulate   <h2> <Mc> <seed>
#   Rscript --vanilla R/gwas/scenarios/resolution_gwas.R resolution <h2> <Mc> <seed> <n> [mean_before|replicates]
#   e.g.  ... simulate 0.6 5 1   then   ... resolution 0.6 5 1 6
# -----------------------------------------------------------------------------
# ARGUMENTS
#   step      "simulate" or "resolution"
#   h2        numeric in [0, 1]   share of the latent genetic variance from the causal SNPs
#   Mc        integer             causal SNPs per parameter
#   seed      integer             seed of the population (seed + 1, seed + 2 for norms, measurement)
#   n         integer             measurement points (step "resolution" only)
#   handling  "mean_before" (default) or "replicates"   replicate handling of the PPIs
# All other settings are the defaults of sim_config() in R/sim/simulation_config.R.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

args <- commandArgs(trailingOnly = TRUE)
usage <- paste("usage: resolution_gwas.R simulate <h2> <Mc> <seed>\n",
               "      resolution_gwas.R resolution <h2> <Mc> <seed> <n> [mean_before|replicates]")
if (length(args) < 4 || !args[1] %in% c("simulate", "resolution")) stop(usage)
step <- args[1]
h2 <- as.numeric(args[2])
Mc <- as.integer(args[3])
seed <- as.integer(args[4])
if (anyNA(c(h2, Mc, seed))) stop(usage)

source(here::here("R", "pipeline", "load_pipeline.R"))
POPULATION_DIR <- here::here("output", "populations", sprintf("resolution_%g_%d_seed%d", round(100 * h2), Mc, seed))

if (step == "simulate") {
  t0 <- Sys.time()
  simulate_population(sim_config(h2 = h2, Mc = Mc, seed = seed), POPULATION_DIR)   # -> data/, variance/
  message("Population simulated in ", format(round(Sys.time() - t0, 1)), ": ", POPULATION_DIR)
} else {
  if (length(args) < 5) stop(usage)
  n_points <- as.integer(args[5])
  replicate_handling <- if (length(args) >= 6) args[6] else "mean_before"
  if (is.na(n_points) || n_points < 2) stop("n must be an integer >= 2")
  if (!file.exists(file.path(POPULATION_DIR, "data", "norms.rds"))) {
    stop("no simulated population in ", POPULATION_DIR, " (run the step simulate first)")
  }
  t0 <- Sys.time()
  RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, n_points)              # -> n_<n>/data/, variance/
  compute_population_ppis(RESOLUTION_DIR, replicate_handling)                 # -> n_<n>/<handling>/ppis/
  message("PPIs done after ", format(round(Sys.time() - t0, 1)))
  run_population_gwas(RESOLUTION_DIR, replicate_handling)                     # -> n_<n>/<handling>/gwas/, plots/
  message("Resolution ", n_points, " done in ", format(round(Sys.time() - t0, 1)), ": ", RESOLUTION_DIR)
}
