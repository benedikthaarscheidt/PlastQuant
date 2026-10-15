# =============================================================================
# compute_population_ppis.R — Step 3: the PPIs of a resolution
# =============================================================================
# Steps on a population folder (see population_folder.R and the README):
#   simulate_population() -> sample_resolution() -> compute_population_ppis() -> run_population_gwas()

# PPIs of a resolution folder (its data/) -> dir/ppis (<PPI>.csv and ppi_scores.rds).
# replicate_handling: "mean_before" (replicates averaged before the PPI) or
# "replicates" (one value per replicate reaction norm; see compute_ppis()).
compute_population_ppis <- function(dir, registry = default_ppi_registry(), exclude = NULL,
                                    replicate_handling = c("mean_before", "replicates")) {
  replicate_handling <- match.arg(replicate_handling)
  pop <- read_population(dir, with_snps = FALSE)
  traits <- compute_ppis(pop$phenotypes, pop$environments, registry, exclude, replicate_handling)
  unlink(file.path(dir, "ppis"), recursive = TRUE)
  ppi_dir <- population_subdir(dir, "ppis")
  saveRDS(traits, file.path(ppi_dir, "ppi_scores.rds"))
  for (n in names(traits)) write.csv(traits[[n]], file.path(ppi_dir, paste0(n, ".csv")), row.names = FALSE)
  message("PPIs written to ", ppi_dir)
  dir
}
