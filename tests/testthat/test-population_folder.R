new_dir <- function(name) {
  d <- file.path(tempdir(), name)
  unlink(d, recursive = TRUE)
  d
}

test_that("a trait matrix with a genotype column becomes the long table and back", {
  tab <- data.frame(Genotype = c("G1", "G1", "G2"), E1 = c(1, 2, 3), E2 = c(4, 5, 6))
  long <- trait_matrix_to_long(tab, c(E1 = 1, E2 = 10))
  expect_equal(names(long), c("Genotype", "Replicate", "Environment", "Trait"))
  expect_equal(long$Replicate[long$Environment == 1], c(1L, 2L, 1L))
  expect_equal(long$Trait[long$Genotype == "G1" & long$Replicate == 2], c(2, 5))
  back <- long_to_trait_matrix(long, c(E1 = 1, E2 = 10))
  expect_equal(back, tab, ignore_attr = TRUE)
})

test_that("a simulated population is written to data/ and read back unchanged", {
  sim <- simulate_reaction_norms(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 6, n_snps = 300))
  dir <- write_population(sim, new_dir("pop_sim"))
  expect_true(all(file.exists(file.path(dir, "data", c("traits/linear.csv", "traits/wave.csv", "environments.csv",
                                                        "snp_matrix.rds", "truth.csv", "config.rds",
                                                        "mother_parameters.csv", "simulation.rds")))))
  pop <- read_population(dir)
  expect_equal(names(pop$phenotypes), sim$config$shapes)
  expect_equal(pop$phenotypes$gaussian, sim$phenotypes$gaussian, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(pop$environments, sim$environments, ignore_attr = TRUE)
  expect_equal(pop$snp_matrix, sim$snp_matrix)
  expect_equal(pop$truth[, c("Shape", "Parameter", "SNP")], sim$truth[, c("Shape", "Parameter", "SNP")], ignore_attr = TRUE)
  expect_equal(pop$config$h2, 0.4)
  expect_equal(nrow(pop$mother_parameters), 24L)
})

test_that("own data: CSV files placed in data/ are read as a population", {
  dir <- new_dir("pop_own")
  dir.create(file.path(dir, "data", "traits"), recursive = TRUE)
  write.csv(data.frame(Genotype = c(101, 101, 102, 103), low = c(1, 1.1, 2, 3), mid = c(2, 2.1, 2.5, 3), high = c(3, 3.2, 3, 3.1)),
            file.path(dir, "data", "traits", "biomass.csv"), row.names = FALSE)
  write.csv(data.frame(Column = c("low", "mid", "high"), Value = c(5, 15, 25)),
            file.path(dir, "data", "environments.csv"), row.names = FALSE)
  X <- matrix(c(0, 1, 2, 2, 1, 0), 3, dimnames = list(NULL, c("s1", "s2")))
  write.csv(data.frame(Genotype = c(101, 102, 103), X, check.names = FALSE), file.path(dir, "data", "snp_matrix.csv"), row.names = FALSE)
  K <- diag(3); dimnames(K) <- list(c("101", "102", "103"), c("101", "102", "103"))
  write.csv(data.frame(Genotype = rownames(K), K, check.names = FALSE), file.path(dir, "data", "kinship.csv"), row.names = FALSE)
  write.csv(data.frame(Trait = "biomass", Parameter = "known", SNP = "s2"), file.path(dir, "data", "truth.csv"), row.names = FALSE)

  pop <- read_population(dir)
  expect_equal(names(pop$phenotypes), "biomass")
  expect_equal(pop$environments, c(5, 15, 25))
  ph <- pop$phenotypes$biomass
  expect_equal(ph$Replicate[ph$Genotype == "101" & ph$Environment == 5], c(1L, 2L))
  expect_equal(rownames(pop$snp_matrix), c("101", "102", "103"))
  expect_equal(pop$kinship, K)
  expect_equal(pop$truth$Shape, "biomass")
  expect_null(pop$config)
})

test_that("own data without SNPs, kinship and truth is allowed", {
  dir <- new_dir("pop_min")
  dir.create(file.path(dir, "data", "traits"), recursive = TRUE)
  write.csv(data.frame(Genotype = c("a", "b"), E1 = 1:2, E2 = 3:4), file.path(dir, "data", "traits", "height.csv"), row.names = FALSE)
  write.csv(data.frame(Column = c("E1", "E2"), Value = c(1, 2)), file.path(dir, "data", "environments.csv"), row.names = FALSE)
  pop <- read_population(dir)
  expect_null(pop$snp_matrix)
  expect_null(pop$kinship)
  expect_null(pop$truth)
})

test_that("format errors name the file and what is wrong", {
  dir <- new_dir("pop_bad")
  expect_error(read_population(dir), "no data/traits")
  dir.create(file.path(dir, "data", "traits"), recursive = TRUE)
  write.csv(data.frame(Genotype = "a", E1 = 1, E9 = 2), file.path(dir, "data", "traits", "height.csv"), row.names = FALSE)
  expect_error(read_population(dir), "environments.csv")
  write.csv(data.frame(Column = "E1", Value = 1), file.path(dir, "data", "environments.csv"), row.names = FALSE)
  expect_error(read_population(dir), "height.csv.*E9")
  write.csv(data.frame(ID = "a", E1 = 1), file.path(dir, "data", "traits", "height.csv"), row.names = FALSE)
  expect_error(read_population(dir), "height.csv.*Genotype")
})

test_that("analyses on dense curves stop when the population has too few points", {
  dir <- write_population(simulate_reaction_norms(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 4, n_snps = 200)),
                          new_dir("pop_points"))
  expect_error(require_population_points(dir, 50), "at least 50 measurement points.*has 6")
  expect_equal(require_population_points(dir, 6), 6L)
  expect_error(require_population_points(new_dir("pop_none"), 6), "no population")
})
