# =============================================================================
# population_steps.R — the steps on a population folder
# =============================================================================
# A population folder holds the population in data/ (see population_folder.R) and one
# subfolder per sampling resolution. Each step reads what it needs and (re)writes only
# its own folder:
#   simulate_population(cfg, dir)     -> data/, variance/      true mother and replicate norms
#   sample_resolution(dir, n)         -> n_<n>/data/, n_<n>/variance/   measurement at n points
#   compute_population_ppis(res_dir)  -> n_<n>/ppis/           plasticity indices
#   run_population_gwas(res_dir)      -> n_<n>/gwas/           GWAS of the PPIs (optional, needs SNPs)
# The figure and analysis scripts read a resolution folder (RESOLUTION_DIR) and write its
# plots/ and analyses/; scripts across resolutions read all n_* folders of POPULATION_DIR.

# Simulate the true reaction norms of a population and save them as a population
# folder: data/ (norms, parameters, SNPs, truth, settings) and variance/ (latent level)
simulate_population <- function(cfg, dir) {
  norms <- simulate_norms(cfg)
  write_population(norms, dir)
  unlink(file.path(dir, "variance"), recursive = TRUE)
  var_dir <- population_subdir(dir, "variance")
  write.csv(latent_variance(norms), file.path(var_dir, "variance_latent.csv"), row.names = FALSE)
  message("Population written to ", dir)
  dir
}

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

# PPIs of a resolution folder (its data/) -> dir/ppis (<PPI>.csv and ppi_scores.rds).
# replicate_handling: "mean_before" (replicates averaged before the PPI) or
# "replicates" (one value per replicate reaction norm; see compute_ppis()).
compute_population_ppis <- function(dir, registry = default_ppi_registry(), exclude = NULL,
                                    replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  pop <- read_population(dir)
  traits <- compute_ppis(pop$phenotypes, pop$environments, registry, exclude, replicate_handling)
  unlink(file.path(dir, "ppis"), recursive = TRUE)
  ppi_dir <- population_subdir(dir, "ppis")
  saveRDS(traits, file.path(ppi_dir, "ppi_scores.rds"))
  for (n in names(traits)) write.csv(traits[[n]], file.path(ppi_dir, paste0(n, ".csv")), row.names = FALSE)
  message("PPIs written to ", ppi_dir)
  dir
}

# GWAS of every PPI in the resolution folder's ppis/ against the population's SNPs ->
# dir/gwas/gwas_results.csv (kinship from the population's data/kinship.csv if present,
# otherwise from the SNPs); with a truth table also the power per shape -> gwas_power.csv
run_population_gwas <- function(dir) {
  ppi_file <- file.path(dir, "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) stop("no PPIs in ", dir, " (run compute_population_ppis(dir) first)")
  pop <- read_population(dir)
  if (is.null(pop$snp_matrix)) stop("the GWAS needs a SNP matrix in ", file.path(pop$population_dir, "data"), " (snp_matrix.csv)")
  traits <- readRDS(ppi_file)
  unlink(file.path(dir, "gwas"), recursive = TRUE)
  gwas_dir <- population_subdir(dir, "gwas")
  kin <- if (is.null(pop$kinship)) kinship_matrix(pop$snp_matrix) else pop$kinship
  res <- run_gwas(traits, pop$snp_matrix, kin)
  write.csv(res, file.path(gwas_dir, "gwas_results.csv"), row.names = FALSE)
  if (!is.null(pop$truth)) {
    power <- evaluate_power(res, pop$truth, ncol(pop$snp_matrix))
    write.csv(power, file.path(gwas_dir, "gwas_power.csv"), row.names = FALSE)
  }
  message("GWAS written to ", gwas_dir)
  dir
}
