truth <- data.frame(Shape = c("linear", "linear", "linear", "gaussian"), Parameter = c("L", "a", "a", "w"),
                    SNP = c("SNP_1", "SNP_2", "SNP_3", "SNP_9"), stringsAsFactors = FALSE)

test_that("power is computed against the causal SNPs of the same shape", {
  hits <- data.frame(Score = "RN", SNP = c("SNP_2", "SNP_9", "SNP_50"), stringsAsFactors = FALSE)
  res <- compute_power(hits, truth[truth$Shape == "linear", ], "RN")
  expect_setequal(res$Parameter, c("L", "a", "False_Positives"))
  expect_equal(res$Power[res$Parameter == "a"], 0.5)
  expect_equal(res$Found[res$Parameter == "False_Positives"], 2)
})

test_that("evaluate_power applies Bonferroni and splits by shape", {
  g <- data.frame(SNP = c("SNP_2", "SNP_9", "SNP_9"), P_Value = c(1e-9, 1e-9, 0.01),
                  Score = "RN", Type = c("linear", "gaussian", "linear"), stringsAsFactors = FALSE)
  res <- evaluate_power(g, truth, n_snps = 1000)
  expect_equal(res$Power[res$Shape == "gaussian" & res$Parameter == "w"], 1)
  expect_equal(res$Found[res$Shape == "linear" & res$Parameter == "False_Positives"], 0)
})
