# =============================================================================
# simulate_population.R — Step 1 of a population folder: simulate the true reaction norms
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# Simulate the true reaction norms of a population and save them as a population
# folder: data/ (norms, parameters, SNPs, truth, settings, plot of the true norms
# reaction_norms.pdf) and variance/ (latent level)
simulate_population <- function(cfg, dir) {
  norms <- simulate_norms(cfg)
  write_population(norms, dir)
  unlink(file.path(dir, "variance"), recursive = TRUE)
  var_dir <- population_subdir(dir, "variance")
  write.csv(latent_variance(norms), file.path(var_dir, "variance_latent.csv"), row.names = FALSE)
  plot_true_norms(norms, file.path(dir, "data", "reaction_norms.pdf"))
  message("Population written to ", dir)
  dir
}
