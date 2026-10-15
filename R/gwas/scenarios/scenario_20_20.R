# =============================================================================
# scenario_20_20.R — heritability scenario: h2 = 0.2, 20 causal SNPs per parameter
# =============================================================================
# Simulates the population (R/sim) into output/populations/scenario_20_20/ (data/,
# variance/), measures it at 6 points (n_6/), computes all PPIs from the replicate-averaged
# reaction norms (n_6/mean_before/ppis/) and runs the GWAS with the power per shape and
# its heatmaps (n_6/mean_before/gwas/, n_6/mean_before/plots/).
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN: Rscript --vanilla R/gwas/scenarios/scenario_20_20.R   (from the repo root)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

source(here::here("R", "pipeline", "load_pipeline.R"))
POPULATION_DIR <- here::here("output", "populations", "scenario_20_20")
simulate_population(sim_config(h2 = 0.2, Mc = 20), POPULATION_DIR)   # true reaction norms -> data/, variance/
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # measured at 6 points -> n_6/data/
compute_population_ppis(RESOLUTION_DIR)                   # PPIs -> n_6/mean_before/ppis/
run_population_gwas(RESOLUTION_DIR)                       # GWAS, power, heatmaps -> n_6/mean_before/gwas/, plots/
