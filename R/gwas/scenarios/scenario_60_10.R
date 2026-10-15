# =============================================================================
# scenario_60_10.R — heritability scenario: h2 = 0.6, 10 causal SNPs per parameter
# =============================================================================
# Simulates the population (R/sim) into output/populations/scenario_60_10/ (data/,
# variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
# GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
# (n_6/plots/, 05_evaluate_gwas.R).
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN: Rscript --vanilla R/gwas/scenarios/scenario_60_10.R   (from the repo root)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

source(here::here("R", "pipeline", "load_pipeline.R"))
POPULATION_DIR <- here::here("output", "populations", "scenario_60_10")
simulate_population(sim_config(h2 = 0.6, Mc = 10), POPULATION_DIR)   # true reaction norms -> data/, variance/
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # measured at 6 points -> n_6/data/
compute_population_ppis(RESOLUTION_DIR)                   # PPIs -> n_6/ppis/
run_population_gwas(RESOLUTION_DIR)                       # GWAS of the PPIs -> n_6/gwas/

# Power heatmaps per shape at this resolution -> n_6/plots/
source(here::here("R", "gwas", "05_evaluate_gwas.R"))
