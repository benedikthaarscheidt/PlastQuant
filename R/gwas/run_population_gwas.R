# =============================================================================
# run_population_gwas.R — Step 4 (optional): the GWAS of the PPIs of a resolution
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# GWAS of every PPI of a resolution and replicate handling (dir/<replicate_handling>/ppis)
# against the population's SNPs -> dir/<replicate_handling>/gwas/gwas_results.csv
# (kinship from the population's data/kinship.csv if present, otherwise from the SNPs);
# with a truth table also the power per shape -> gwas_power.csv
run_population_gwas <- function(dir, replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  set_dir <- ppi_set_path(dir, replicate_handling)
  ppi_file <- file.path(set_dir, "ppis", "ppi_scores.rds")
  if (!file.exists(ppi_file)) {
    stop("no PPIs in ", set_dir, " (run compute_population_ppis(dir, replicate_handling = \"", replicate_handling, "\") first)")
  }
  pop <- read_population(dir)
  if (is.null(pop$snp_matrix)) stop("the GWAS needs a SNP matrix in ", file.path(pop$population_dir, "data"), " (snp_matrix.csv)")
  traits <- readRDS(ppi_file)
  gwas_dir <- file.path(set_dir, "gwas")
  unlink(gwas_dir, recursive = TRUE)
  dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
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
