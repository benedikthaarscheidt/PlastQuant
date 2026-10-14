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

test_that("the population holds the true norms; a resolution folder holds one measurement of them", {
  norms <- simulate_norms(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 6, n_snps = 300))
  dir <- write_population(norms, new_dir("pop_sim"))
  expect_true(all(file.exists(file.path(dir, "data", c("norms.rds", "mother_parameters.csv", "plant_parameters.csv",
                                                        "snp_matrix.rds", "truth.csv", "config.rds")))))
  expect_false(dir.exists(file.path(dir, "data", "traits")))
  expect_error(read_population(dir), "not a resolution folder")

  measured <- measure_norms(norms, 6)
  res <- write_resolution(measured, resolution_path(dir, 6))
  expect_equal(basename(res), "n_6")
  expect_true(all(file.exists(file.path(res, "data", c("traits/linear.csv", "traits/wave.csv", "environments.csv")))))
  pop <- read_population(res)
  expect_equal(names(pop$phenotypes), norms$config$shapes)
  expect_equal(pop$phenotypes$gaussian, measured$phenotypes$gaussian, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(pop$environments, measured$environments, ignore_attr = TRUE)
  expect_equal(pop$snp_matrix, norms$snp_matrix)                # from the population's data/
  expect_equal(pop$truth$SNP, norms$truth$SNP)
  expect_equal(pop$config$h2, 0.4)
  expect_equal(nrow(pop$mother_parameters), 24L)
  expect_equal(normalizePath(pop$population_dir), normalizePath(dir))
})

test_that("own data: SNPs, kinship and truth in data/, the measurements in a resolution folder", {
  dir <- new_dir("pop_own")
  dir.create(file.path(dir, "data"), recursive = TRUE)
  dir.create(file.path(dir, "n_3", "data", "traits"), recursive = TRUE)
  write.csv(data.frame(Genotype = c(101, 101, 102, 103), low = c(1, 1.1, 2, 3), mid = c(2, 2.1, 2.5, 3), high = c(3, 3.2, 3, 3.1)),
            file.path(dir, "n_3", "data", "traits", "biomass.csv"), row.names = FALSE)
  write.csv(data.frame(Column = c("low", "mid", "high"), Value = c(5, 15, 25)),
            file.path(dir, "n_3", "data", "environments.csv"), row.names = FALSE)
  X <- matrix(c(0, 1, 2, 2, 1, 0), 3, dimnames = list(NULL, c("s1", "s2")))
  write.csv(data.frame(Genotype = c(101, 102, 103), X, check.names = FALSE), file.path(dir, "data", "snp_matrix.csv"), row.names = FALSE)
  K <- diag(3); dimnames(K) <- list(c("101", "102", "103"), c("101", "102", "103"))
  write.csv(data.frame(Genotype = rownames(K), K, check.names = FALSE), file.path(dir, "data", "kinship.csv"), row.names = FALSE)
  write.csv(data.frame(Trait = "biomass", Parameter = "known", SNP = "s2"), file.path(dir, "data", "truth.csv"), row.names = FALSE)

  pop <- read_population(file.path(dir, "n_3"))
  expect_equal(names(pop$phenotypes), "biomass")
  expect_equal(pop$environments, c(5, 15, 25))
  ph <- pop$phenotypes$biomass
  expect_equal(ph$Replicate[ph$Genotype == "101" & ph$Environment == 5], c(1L, 2L))
  expect_equal(rownames(pop$snp_matrix), c("101", "102", "103"))
  expect_equal(pop$kinship, K)
  expect_equal(pop$truth$Shape, "biomass")
  expect_null(pop$config)
  expect_equal(basename(resolution_dirs(dir)), "n_3")
})

test_that("own data without SNPs, kinship and truth is allowed", {
  dir <- new_dir("pop_min")
  dir.create(file.path(dir, "n_2", "data", "traits"), recursive = TRUE)
  write.csv(data.frame(Genotype = c("a", "b"), E1 = 1:2, E2 = 3:4), file.path(dir, "n_2", "data", "traits", "height.csv"), row.names = FALSE)
  write.csv(data.frame(Column = c("E1", "E2"), Value = c(1, 2)), file.path(dir, "n_2", "data", "environments.csv"), row.names = FALSE)
  pop <- read_population(file.path(dir, "n_2"))
  expect_null(pop$snp_matrix)
  expect_null(pop$kinship)
  expect_null(pop$truth)
})

test_that("format errors name the file and what is wrong", {
  res <- file.path(new_dir("pop_bad"), "n_2")
  expect_error(read_population(res), "not a resolution folder")
  dir.create(file.path(res, "data", "traits"), recursive = TRUE)
  expect_error(read_population(res), "no data/traits")
  write.csv(data.frame(Genotype = "a", E1 = 1, E9 = 2), file.path(res, "data", "traits", "height.csv"), row.names = FALSE)
  expect_error(read_population(res), "environments.csv")
  write.csv(data.frame(Column = "E1", Value = 1), file.path(res, "data", "environments.csv"), row.names = FALSE)
  expect_error(read_population(res), "height.csv.*E9")
  write.csv(data.frame(ID = "a", E1 = 1), file.path(res, "data", "traits", "height.csv"), row.names = FALSE)
  expect_error(read_population(res), "height.csv.*Genotype")
})

test_that("the saved norms of a simulated population are sampled at any resolution", {
  cfg <- sim_config(h2 = 0.4, Mc = 2, n_per_shape = 5, n_snps = 200, seed = 7)
  norms <- simulate_norms(cfg)
  dir <- write_population(norms, new_dir("pop_resample"))
  expect_equal(sample_population(dir, 6)$phenotypes, measure_norms(norms, 6)$phenotypes)
  s3 <- sample_population(dir, 3)
  expect_equal(s3$environments, c(1, 5.5, 10))
  expect_equal(unique(s3$phenotypes$gaussian$Genotype), unique(norms$plant_parameters$gaussian$Genotype))
  expect_error(sample_population(new_dir("pop_own_nosim"), 3), "simulated population")
})

test_that("resolution folders of a population, ordered by their number of points", {
  norms <- simulate_norms(sim_config(h2 = 0.4, Mc = 2, n_per_shape = 3, n_snps = 200))
  dir <- write_population(norms, new_dir("pop_res_list"))
  for (n in c(11, 3, 6)) write_resolution(measure_norms(norms, n), resolution_path(dir, n))
  dir.create(file.path(dir, "analyses"))
  expect_equal(basename(resolution_dirs(dir)), c("n_3", "n_6", "n_11"))
})

test_that("own data built in R is written as population plus one resolution folder", {
  ph <- data.frame(Genotype = rep(c("A", "B"), each = 3), Replicate = 1L, Environment = rep(c(1, 2, 3), 2), Trait = 1:6)
  X <- matrix(c(0, 2, 1, 1), 2, dimnames = list(c("A", "B"), c("s1", "s2")))
  dir <- write_population(pipeline_input(list(biomass = ph), c(1, 2, 3), X), new_dir("pop_own_r"))
  expect_true(file.exists(file.path(dir, "data", "snp_matrix.csv")))
  expect_equal(basename(resolution_dirs(dir)), "n_3")
  expect_equal(read_population(resolution_dirs(dir))$phenotypes$biomass$Trait, 1:6, ignore_attr = TRUE)
})
