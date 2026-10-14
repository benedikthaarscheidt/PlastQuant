# =============================================================================
# scenario_40_10.R — heritability scenario: h2 = 0.4, 10 causal SNPs per parameter
# =============================================================================
# Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
# per-shape evaluation, and writes everything to output/scenario_40_10/run_<time>/.
# All other settings (variance components, study space, ...) are the defaults of
# sim_config() in R/sim/simulation_config.R.
# HOW TO RUN: Rscript --vanilla R/scenarios/scenario_40_10.R   (from the repo root)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

source(here::here("R", "pipeline", "load_pipeline.R"))
run_simulation_scenario(sim_config(h2 = 0.4, Mc = 10), out_dir = here::here("output", "scenario_40_10"))
