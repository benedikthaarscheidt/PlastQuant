# =============================================================================
# run_pipeline.R — PPIs -> GWAS -> evaluation on a population folder
# =============================================================================
# A population folder holds the population in data/ (see population_folder.R) and
# everything computed from it; each step rewrites its own subfolder.

# Simulate a population and save it as a population folder: data/ and variance/
simulate_population <- function(cfg, dir) {
  sim <- simulate_reaction_norms(cfg)
  write_population(sim, dir)
  vr <- variance_report(sim)
  unlink(file.path(dir, "variance"), recursive = TRUE)
  var_dir <- population_subdir(dir, "variance")
  write.csv(vr$latent, file.path(var_dir, "variance_latent.csv"), row.names = FALSE)
  write.csv(vr$trait, file.path(var_dir, "variance_trait.csv"), row.names = FALSE)
  write.csv(vr$measurement, file.path(var_dir, "measurement.csv"), row.names = FALSE)
  message("Population written to ", dir)
  dir
}

# PPIs of the population in dir/data -> dir/ppis; with SNPs the GWAS -> dir/gwas
# (kinship from data/kinship.csv if present, otherwise from the SNPs); with a truth
# table also the power per shape -> dir/gwas/gwas_power.csv. Returns dir.
# gwas = NULL runs the GWAS whenever the population has a SNP matrix.
run_pipeline <- function(dir, registry = default_ppi_registry(), exclude = NULL, gwas = NULL,
                         replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  pop <- read_population(dir)
  if (is.null(gwas)) gwas <- !is.null(pop$snp_matrix)
  if (gwas && is.null(pop$snp_matrix)) stop("gwas = TRUE needs a SNP matrix in ", file.path(dir, "data"))

  traits <- compute_ppis(pop$phenotypes, pop$environments, registry, exclude, replicate_handling)
  unlink(file.path(dir, "ppis"), recursive = TRUE)
  ppi_dir <- population_subdir(dir, "ppis")
  saveRDS(traits, file.path(ppi_dir, "ppi_scores.rds"))
  for (n in names(traits)) write.csv(traits[[n]], file.path(ppi_dir, paste0(n, ".csv")), row.names = FALSE)

  unlink(file.path(dir, "gwas"), recursive = TRUE)
  if (gwas) {
    gwas_dir <- population_subdir(dir, "gwas")
    kin <- if (is.null(pop$kinship)) kinship_matrix(pop$snp_matrix) else pop$kinship
    res <- run_gwas(traits, pop$snp_matrix, kin)
    write.csv(res, file.path(gwas_dir, "gwas_results.csv"), row.names = FALSE)
    if (!is.null(pop$truth)) {
      power <- evaluate_power(res, pop$truth, ncol(pop$snp_matrix))
      write.csv(power, file.path(gwas_dir, "gwas_power.csv"), row.names = FALSE)
    }
  }
  message("Results written to ", dir)
  dir
}

# One heritability scenario: simulate the population into dir, then run the pipeline on it
run_simulation_scenario <- function(cfg, dir, registry = default_ppi_registry(), exclude = NULL, gwas = TRUE,
                                    replicate_handling = "mean_before") {
  simulate_population(cfg, dir)
  run_pipeline(dir, registry, exclude, gwas, replicate_handling)
}
