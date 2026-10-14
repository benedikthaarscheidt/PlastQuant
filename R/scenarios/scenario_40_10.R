# =============================================================================
# scenario_40_10.R — heritability scenario: h2 = 0.4, 10 causal SNPs per parameter
# =============================================================================
# Simulates the population (R/sim) and saves it as the population folder
# output/populations/scenario_40_10/ (data/, variance/), computes all PPIs (ppis/),
# runs the GWAS and the per-shape evaluation (gwas/) and draws the power heatmaps
# (plots/, 05_evaluate_gwas.R) into the same folder.
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN: Rscript --vanilla R/scenarios/scenario_40_10.R   (from the repo root)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

source(here::here("R", "pipeline", "load_pipeline.R"))
POPULATION_DIR <- run_simulation_scenario(sim_config(h2 = 0.4, Mc = 10), here::here("output", "populations", "scenario_40_10"))

# Power heatmaps per shape of this population
source(here::here("R", "05_evaluate_gwas.R"))
