# =============================================================================
# run_pipeline.R — PPIs -> GWAS -> evaluation for simulated or own data
# =============================================================================

# input: pipeline_input() or simulate_reaction_norms(); writes everything into a
# new folder out_dir/run_<timestamp> and returns its path
run_pipeline <- function(input, out_dir, registry = default_ppi_registry(), exclude = NULL, gwas = TRUE,
                         replicate_handling = "mean_before") {
  check_pipeline_input(input)
  if (gwas && is.null(input$snp_matrix)) stop("gwas = TRUE needs a snp_matrix; use gwas = FALSE for phenotypes only")
  run_dir <- file.path(out_dir, paste0("run_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
  dir.create(file.path(run_dir, "ppi_scores"), recursive = TRUE)

  traits <- compute_ppis(input$phenotypes, input$environments, registry, exclude, replicate_handling)
  saveRDS(input$phenotypes, file.path(run_dir, "phenotypes.rds"))
  if (!is.null(input$snp_matrix)) saveRDS(input$snp_matrix, file.path(run_dir, "snp_matrix.rds"))
  saveRDS(traits, file.path(run_dir, "ppi_scores.rds"))
  for (n in names(traits)) write.csv(traits[[n]], file.path(run_dir, "ppi_scores", paste0(n, ".csv")), row.names = FALSE)
  if (!is.null(input$truth)) write.csv(input$truth, file.path(run_dir, "truth.csv"), row.names = FALSE)

  if (gwas) {
    kin <- kinship_matrix(input$snp_matrix)
    res <- run_gwas(traits, input$snp_matrix, kin)
    write.csv(res, file.path(run_dir, "gwas_results.csv"), row.names = FALSE)
    if (!is.null(input$truth)) {
      power <- evaluate_power(res, input$truth, ncol(input$snp_matrix))
      write.csv(power, file.path(run_dir, "gwas_power.csv"), row.names = FALSE)
    }
  }
  message("Results written to ", run_dir)
  run_dir
}

# One heritability scenario: simulate, report the variance components, run the pipeline
run_simulation_scenario <- function(cfg, out_dir, registry = default_ppi_registry(), exclude = NULL, gwas = TRUE,
                                    replicate_handling = "mean_before") {
  sim <- simulate_reaction_norms(cfg)
  run_dir <- run_pipeline(sim, out_dir, registry, exclude, gwas, replicate_handling)
  vr <- variance_report(sim)
  write.csv(vr$latent, file.path(run_dir, "variance_latent.csv"), row.names = FALSE)
  write.csv(vr$trait, file.path(run_dir, "variance_trait.csv"), row.names = FALSE)
  write.csv(vr$measurement, file.path(run_dir, "measurement.csv"), row.names = FALSE)
  write.csv(sim$mother_parameters, file.path(run_dir, "mother_parameters.csv"), row.names = FALSE)
  saveRDS(sim$mother_curves, file.path(run_dir, "mother_curves.rds"))
  saveRDS(cfg, file.path(run_dir, "config.rds"))
  run_dir
}
