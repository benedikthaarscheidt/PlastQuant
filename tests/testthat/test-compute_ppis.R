toy <- function() {
  env <- c(1, 4, 7, 10)
  d <- expand.grid(Environment = env, Replicate = 1:2, Genotype = paste0("G", 1:3), stringsAsFactors = FALSE)
  d$Trait <- 1 + as.integer(sub("G", "", d$Genotype)) * d$Environment / 10 + (d$Replicate - 1.5) * 0.01
  list(phenotypes = list(shapeA = d[, c("Genotype", "Replicate", "Environment", "Trait")]), env = env)
}
slope <- function(trait_values, env_values) unname(coef(lm(trait_values ~ env_values))[2])
reg <- list(
  slope = ppi_entry("slope", slope, inputs = c(trait_values = "trait", env_values = "env")),
  mean  = ppi_entry("mean", function(trait_values) mean(trait_values)),
  pop   = ppi_entry("pop", function(Y) rowMeans(Y), inputs = c(Y = "matrix"), level = "population",
                    extract = function(res) res)
)

test_that("replicate means form a genotypes x environments matrix", {
  x <- toy()
  M <- replicate_means(x$phenotypes$shapeA, x$env)
  expect_equal(dim(M), c(3L, 4L))
  expect_equal(rownames(M), c("G1", "G2", "G3"))
})

test_that("each PPI gets its declared inputs; result layout matches scores_list", {
  x <- toy()
  sc <- compute_ppis(x$phenotypes, x$env, reg)
  expect_named(sc, c("slope", "mean", "pop"))
  expect_equal(names(sc$slope), c("Score", "Type", "Genotype_ID", "Replicate"))
  expect_equal(sc$slope$Score, c(0.1, 0.2, 0.3), tolerance = 1e-10)
  expect_equal(sc$pop$Score, sc$mean$Score)
  expect_equal(unique(sc$mean$Type), "shapeA")
})

test_that("PPIs can be excluded; unknown names are an error", {
  x <- toy()
  expect_named(compute_ppis(x$phenotypes, x$env, reg, exclude = c("pop", "mean")), "slope")
  expect_error(compute_ppis(x$phenotypes, x$env, reg, exclude = "nope"), "cannot exclude unknown PPI")
})

test_that("a failing PPI gives NA instead of stopping the run", {
  x <- toy()
  bad <- list(bad = ppi_entry("bad", function(trait_values) stop("boom")))
  expect_message(sc <- compute_ppis(x$phenotypes, x$env, bad), "no value for 3")
  expect_true(all(is.na(sc$bad$Score)))
})

test_that("the default registry runs on simulated data", {
  sim <- simulate_reaction_norms(sim_config(h2 = 0.4, Mc = 5, n_per_shape = 8, n_snps = 500, seed = 2))
  sc <- suppressMessages(compute_ppis(sim$phenotypes, sim$environments))
  expect_equal(length(sc), 28L)
  expect_equal(nrow(sc$FW), 32L)
})

test_that("replicate handling: averaged before, or one row per replicate", {
  x <- toy()
  reps <- compute_ppis(x$phenotypes, x$env, reg["slope"], replicate_handling = "replicates")$slope
  expect_equal(nrow(reps), 6L)
  expect_equal(reps$Genotype_ID, rep(c("G1", "G2", "G3"), each = 2))
  expect_equal(reps$Replicate, rep(1:2, 3))
  expect_error(compute_ppis(x$phenotypes, x$env, reg, replicate_handling = "mean_after"), "should be one of")
})

test_that("a population-level PPI that errors gives NA with a message", {
  x <- toy()
  bad <- list(fw = ppi_entry("fw", function(Y) stop("singular"), inputs = c(Y = "matrix"), level = "population",
                             extract = function(res) setNames(res$beta, res$genotype)))
  expect_message(sc <- compute_ppis(x$phenotypes, x$env, bad), "no value for 3 genotype")
  expect_true(all(is.na(sc$fw$Score)))
  expect_equal(sc$fw$Genotype_ID, c("G1", "G2", "G3"))
})

test_that("factor and numeric genotype IDs keep their scores", {
  x <- toy()
  d <- x$phenotypes$shapeA
  d <- d[order(match(d$Genotype, c("G3", "G1", "G2"))), ]
  d$Genotype <- factor(d$Genotype)
  sc <- compute_ppis(list(shapeA = d), x$env, reg["slope"], replicate_handling = "replicates")
  expect_equal(sc$slope$Score[sc$slope$Genotype_ID == "G3"], c(0.3, 0.3), tolerance = 1e-10)
  expect_equal(sc$slope$Score[sc$slope$Genotype_ID == "G1"], c(0.1, 0.1), tolerance = 1e-10)
  n <- x$phenotypes$shapeA
  n$Genotype <- 100 + as.integer(sub("G", "", n$Genotype))
  sc <- compute_ppis(list(shapeA = n), x$env, reg["slope"])
  expect_equal(sc$slope$Genotype_ID, c("101", "102", "103"))
  expect_equal(sc$slope$Score, c(0.1, 0.2, 0.3), tolerance = 1e-10)
})

test_that("a population PPI is always computed once on the replicate means, one value per genotype", {
  x <- toy()
  d <- x$phenotypes$shapeA
  d$Trait[d$Replicate == 2] <- d$Trait[d$Replicate == 2] + c(0, 3, 9)[as.integer(sub("G", "", d$Genotype[d$Replicate == 2]))]
  # the result depends on the whole population: centre each row on the population mean
  pop_c <- ppi_entry("pop_c", function(Y) rowMeans(Y) - mean(Y), inputs = c(Y = "matrix"), level = "population",
                     extract = function(res) res)
  before <- compute_ppis(list(shapeA = d), x$env, list(pop_c = pop_c))$pop_c
  reps <- compute_ppis(list(shapeA = d), x$env, list(pop_c = pop_c), replicate_handling = "replicates")$pop_c
  expect_equal(nrow(reps), 3L)
  expect_equal(reps$Genotype_ID, c("G1", "G2", "G3"))
  expect_true(all(is.na(reps$Replicate)))
  expect_equal(reps$Score, before$Score)
})

test_that("PPI results as a wide table: one row per genotype, one column per PPI, rows by shape", {
  traits <- list(
    RN = data.frame(Score = c(1, 2, 3), Type = c("wave", "linear", "linear"), Genotype_ID = c("G9", "G1", "G2"), Replicate = NA),
    CV = data.frame(Score = c(30, 10, 20), Type = c("wave", "linear", "linear"), Genotype_ID = c("G9", "G1", "G2"), Replicate = NA)
  )
  w <- ppi_wide_table(traits, shapes = c("linear", "wave"))
  expect_equal(names(w), c("Shape", "Genotype", "RN", "CV"))
  expect_equal(w$Genotype, c("G1", "G2", "G9"))
  expect_equal(w$CV, c(10, 20, 30))
  expect_error(ppi_wide_table(lapply(traits, function(d) transform(d, Replicate = 1L))), "one value per genotype")
})
