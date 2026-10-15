fake_population <- function(root, scenario, h2, Mc, gwas, truth) {
  dir <- file.path(root, scenario)
  dir.create(file.path(dir, "data"), recursive = TRUE)
  dir.create(file.path(dir, "n_6", "mean_before", "gwas"), recursive = TRUE)
  saveRDS(list(h2 = h2, Mc = Mc), file.path(dir, "data", "config.rds"))
  write.csv(gwas, file.path(dir, "n_6", "mean_before", "gwas", "gwas_results.csv"), row.names = FALSE)
  names(truth)[names(truth) == "Shape"] <- "Trait"
  write.csv(truth, file.path(dir, "data", "truth.csv"), row.names = FALSE)
  dir
}

toy_gwas <- function(p_s1) {
  data.frame(SNP = rep(c("s1", "s2", "s3"), 4), LogP = NA, P_Value = c(p_s1, 0.5, 1e-5, rep(0.5, 9)),
             Score = rep(c("RN", "FW"), each = 6), Type = rep(rep(c("linear", "wave"), each = 3), 2))
}
toy_truth <- data.frame(Shape = c("linear", "wave"), Parameter = "a", SNP = c("s1", "s3"), Effect = 1)

test_that("the populations of the scenarios are combined with their h2 and Mc", {
  root <- file.path(tempdir(), "comb1"); unlink(root, recursive = TRUE)
  fake_population(root, "scenario_20_5", 0.2, 5, toy_gwas(1e-6), toy_truth)
  fake_population(root, "scenario_60_5", 0.6, 5, toy_gwas(0.5), toy_truth)
  comb <- combine_scenario_runs(file.path(root, c("scenario_20_5", "scenario_60_5")), n_points = 6, p_threshold = 1e-3)
  expect_equal(sort(unique(comb$hits$H2)), c(0.2, 0.6))
  expect_true("s1" %in% comb$hits$SNP[comb$hits$H2 == 0.2])
  expect_false("s1" %in% comb$hits$SNP[comb$hits$H2 == 0.6])
  expect_true(all(comb$hits$P_Value < 1e-3))
  expect_equal(nrow(comb$truth), 4L)
  expect_true(all(c("Shape", "Parameter", "SNP") %in% names(comb$truth)))
  expect_equal(nrow(comb$tested), 8L)                          # 2 scenarios x 2 PPIs x 2 shapes
})

test_that("detections are counted per shape: a parameter of one shape never counts for another", {
  root <- file.path(tempdir(), "comb2"); unlink(root, recursive = TRUE)
  fake_population(root, "scenario_40_5", 0.4, 5, toy_gwas(1e-6), toy_truth)
  comb <- combine_scenario_runs(file.path(root, "scenario_40_5"), n_points = 6, p_threshold = 1e-3)
  det <- detection_counts(comb)
  lin <- det[det$Shape == "linear" & det$Score == "RN", ]
  expect_equal(lin$Found[lin$Parameter == "a"], 1L)
  # s3 is causal for the wave shape only: in the linear GWAS it is a false positive
  expect_equal(lin$Found[lin$Parameter == "False_Positives"], 1L)
  wav <- det[det$Shape == "wave" & det$Score == "RN", ]
  expect_equal(wav$Found[wav$Parameter == "a"], 0L)
  expect_true(all(c("H2", "Mc", "Shape", "Parameter", "Score", "True_SNPs", "Found", "Power") %in% names(det)))
  # PPIs without any hit still appear, with power 0
  expect_equal(det$Power[det$Shape == "linear" & det$Score == "FW" & det$Parameter == "a"], 0)
})

test_that("a scenario without GWAS results at that resolution is an error naming the folder", {
  root <- file.path(tempdir(), "comb3"); unlink(root, recursive = TRUE)
  dir.create(file.path(root, "scenario_20_10", "data"), recursive = TRUE)
  expect_error(combine_scenario_runs(file.path(root, "scenario_20_10"), n_points = 6), "scenario_20_10/n_6/mean_before")
})
