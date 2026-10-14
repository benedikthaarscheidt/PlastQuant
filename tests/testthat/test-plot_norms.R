test_that("true norm curves come from the parameters: the replicate curves are the noise-free measurement", {
  norms <- simulate_norms(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200))
  ids <- c("Genotype_6", "Genotype_8")
  cur <- true_norm_curves(norms, "gaussian", ids, n_grid = 6)
  expect_equal(sort(unique(cur$Genotype)), ids)
  expect_equal(sort(unique(cur$Replicate)), 0:3)                 # 0 = mother norm
  noise_free <- measure_norms(norms, 6, error = FALSE)$phenotypes$gaussian
  plants <- cur[cur$Replicate > 0, ]
  expect_equal(plants$Trait, noise_free$Trait[match(paste(plants$Genotype, plants$Replicate, plants$Environment),
                                                   paste(noise_free$Genotype, noise_free$Replicate, noise_free$Environment))])
})

test_that("simulating a population draws its true norms into data/", {
  dir <- file.path(tempdir(), "plot_pop"); unlink(dir, recursive = TRUE)
  simulate_population(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200), dir)
  expect_true(file.exists(file.path(dir, "data", "reaction_norms.pdf")))
})
