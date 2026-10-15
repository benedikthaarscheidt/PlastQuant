# =============================================================================
# sample_resolution.R — Step 2: measure the population at a resolution
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# Measure the population's norms at n_points equidistant points: the resolution folder
# dir/n_<n_points> with data/ and variance/ (trait level, measurement error). An earlier
# folder of this resolution is replaced, including its results. Returns its path.
sample_resolution <- function(dir, n_points) {
  measured <- sample_population(dir, n_points)
  res_dir <- resolution_path(dir, n_points)
  unlink(res_dir, recursive = TRUE)
  write_resolution(measured, res_dir)
  tv <- trait_variance(measured)
  var_dir <- population_subdir(res_dir, "variance")
  write.csv(tv$trait, file.path(var_dir, "variance_trait.csv"), row.names = FALSE)
  write.csv(tv$measurement, file.path(var_dir, "measurement.csv"), row.names = FALSE)
  message("Resolution written to ", res_dir)
  res_dir
}
