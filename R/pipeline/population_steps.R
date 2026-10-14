# =============================================================================
# population_steps.R — the steps on a population folder
# =============================================================================
# A population folder holds the population in data/ (see population_folder.R) and
# everything computed from it. Each step reads what it needs from the folder and
# (re)writes only its own subfolder:
#   simulate_population(cfg, dir)    -> data/, variance/   reaction norms and their sampling
#   compute_population_ppis(dir)     -> ppis/              plasticity indices
#   run_population_gwas(dir)         -> gwas/              GWAS of the PPIs (optional, needs SNPs)
# The figure and analysis scripts read data/ and ppis/ and write plots/ and analyses/.

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

# PPIs of the population in dir/data -> dir/ppis (<PPI>.csv and ppi_scores.rds).
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

# GWAS of every PPI in dir/ppis against the SNPs in dir/data -> dir/gwas/gwas_results.csv
# (kinship from data/kinship.csv if present, otherwise from the SNPs); with a truth
# table also the power per shape -> dir/gwas/gwas_power.csv
run_population_gwas <- function(dir) {
  ppi_file <- file.path(dir, "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) stop("no PPIs in ", dir, " (run compute_population_ppis(dir) first)")
  pop <- read_population(dir)
  if (is.null(pop$snp_matrix)) stop("the GWAS needs a SNP matrix in ", file.path(dir, "data"), " (snp_matrix.csv)")
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
