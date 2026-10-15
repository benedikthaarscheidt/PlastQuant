# =============================================================================
# run_population_gwas.R — Step 4 (optional): the GWAS of the PPIs of a resolution
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# GWAS of every PPI of a resolution and replicate handling (dir/<replicate_handling>/ppis)
# against the population's SNPs -> dir/<replicate_handling>/gwas/gwas_results.csv
# (kinship from the population's data/kinship.csv if present, otherwise from the SNPs);
# with a truth table also the power per shape (Bonferroni alpha / number of SNPs) ->
# gwas_power.csv and its heatmaps -> dir/<replicate_handling>/plots/gwas_power_heatmap.pdf
run_population_gwas <- function(dir, replicate_handling = c("mean_before", "replicates"), alpha = 0.05) {
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
    n_snps <- ncol(pop$snp_matrix)
    power <- evaluate_power(res, pop$truth, n_snps, alpha = alpha)
    write.csv(power, file.path(gwas_dir, "gwas_power.csv"), row.names = FALSE)
    cfg <- pop$config
    subtitle <- paste0(
      if (!is.null(cfg)) paste0("h2 = ", cfg$h2, ", ", cfg$Mc, " causal SNPs per parameter, ",
                                cfg$n_per_shape, " genotypes per shape, "),
      n_snps, " SNPs, ", replicate_handling, "  (Bonferroni P < ", format(alpha / n_snps, digits = 2, scientific = TRUE), ")"
    )
    unlink(file.path(set_dir, "plots", "gwas_power_heatmap.pdf"))
    plot_gwas_power(power, file.path(set_dir, "plots", "gwas_power_heatmap.pdf"), subtitle)
  }
  message("GWAS written to ", gwas_dir)
  dir
}
