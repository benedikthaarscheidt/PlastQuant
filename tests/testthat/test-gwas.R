test_that("GWAS finds a strong causal SNP", {
  set.seed(4)
  X <- simulate_snp_genotypes(200, 300)$X
  y <- setNames(as.numeric(X[, "SNP_7"]) + rnorm(200, 0, 0.3), rownames(X))
  traits <- list(toy = data.frame(Score = unname(y), Type = "linear", Genotype_ID = names(y)))
  res <- suppressMessages(run_gwas(traits, X))
  expect_equal(res$SNP[which.min(res$P_Value)], "SNP_7")
  expect_equal(nrow(res), 300L)
})

test_that("GWAS on replicate values: each line occurs several times", {
  set.seed(5)
  X <- simulate_snp_genotypes(150, 200)$X
  line_value <- 0.8 * X[, "SNP_3"]
  ids <- rep(rownames(X), each = 3)
  y <- line_value[ids] + rnorm(450, 0, 0.5)
  traits <- list(toy = data.frame(Score = unname(y), Type = "linear", Genotype_ID = ids))
  res <- suppressMessages(run_gwas(traits, X))
  expect_equal(res$SNP[which.min(res$P_Value)], "SNP_3")
})
