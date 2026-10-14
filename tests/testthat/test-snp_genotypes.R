test_that("SNP matrix has the right shape, names and values", {
  set.seed(1)
  g <- simulate_snp_genotypes(40, 300)
  expect_equal(dim(g$X), c(40L, 300L))
  expect_equal(rownames(g$X)[c(1, 40)], c("Genotype_1", "Genotype_40"))
  expect_equal(colnames(g$X)[300], "SNP_300")
  expect_true(is.integer(g$X))
  expect_true(all(g$X %in% 0:2))
  expect_true(all(g$p >= 0.2 & g$p <= 0.8))
})

test_that("causal SNPs are disjoint, Mc per parameter incl. V, effects +-1", {
  set.seed(2)
  tr <- draw_causal_snps(sim_config(h2 = 0.4, Mc = 10, n_snps = 5000))
  expect_equal(nrow(tr), 190L)
  expect_equal(anyDuplicated(tr$SNP_index), 0L)
  expect_true(all(table(paste(tr$Shape, tr$Parameter)) == 10))
  expect_equal(tr$Parameter[tr$Shape == "wave"][seq(1, 60, by = 10)], c("L", "a", "k", "c", "w", "V"))
  expect_true(all(tr$Effect %in% c(-1L, 1L)))
})

test_that("too few SNPs is a clear error", {
  expect_error(draw_causal_snps(sim_config(h2 = 0.4, Mc = 20, n_snps = 100)), "need 380 causal SNPs but only 100")
})

test_that("SNP score has mean 0 and variance 1", {
  set.seed(3)
  g <- simulate_snp_genotypes(20000, 50)
  s <- snp_score(g$X, g$p, rows = seq_len(20000), snp_idx = 1:20, beta = rep(c(1L, -1L), 10))
  expect_equal(mean(s), 0, tolerance = 0.03)
  expect_equal(var(s), 1, tolerance = 0.05)
})
