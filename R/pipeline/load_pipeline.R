# =============================================================================
# load_pipeline.R — load all functions of the repository (no script is run)
# =============================================================================
# The folders also hold scripts, so the function files are listed explicitly:
# sourcing this file only defines functions.
pipeline_function_files <- c(
  # simulation model and simulating a population
  "sim/simulation_config.R", "sim/snp_genotypes.R", "sim/genetic_values.R", "sim/parameter_mapping.R",
  "sim/reaction_norm_shapes.R", "sim/measurement.R", "sim/simulate_reaction_norms.R", "sim/variance_report.R",
  "sim/plot_norms.R", "sim/simulate_population.R", "sim/snp_panel.R",
  # infrastructure: data contract, population folders, sampling a resolution
  "pipeline/pipeline_input.R", "pipeline/population_folder.R", "pipeline/sample_resolution.R",
  # PPI computation
  "ppi/ppi_registry.R", "ppi/compute_ppis.R", "ppi/compute_population_ppis.R",
  # GWAS and its evaluation
  "gwas/gwas.R", "gwas/evaluate_power.R", "gwas/plot_gwas_power.R", "gwas/run_population_gwas.R", "gwas/combine_scenario_runs.R",
  # comparison of the PPIs
  "analysis/similarity.R", "analysis/resolution_stability.R"
)
for (f in pipeline_function_files) source(here::here("R", f))
rm(f)
