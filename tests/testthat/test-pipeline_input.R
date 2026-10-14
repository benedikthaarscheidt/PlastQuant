own <- function() {
  X <- matrix(c(0L, 1L, 2L, 1L, 0L, 2L), 3, 2, dimnames = list(c("var_A", "var_B", "var_C"), c("s1", "s2")))
  ph <- data.frame(Genotype = rep(c("var_A", "var_B", "var_C"), each = 2), Replicate = 1L,
                   Environment = rep(c(5, 9), 3), Trait = c(1, 2, 1, 3, 2, 2))
  list(X = X, ph = list(biomass = ph))
}

test_that("own data passes the check", {
  o <- own()
  x <- pipeline_input(o$ph, environments = c(5, 9), snp_matrix = o$X, truth = data.frame(Shape = "biomass", Parameter = "known", SNP = "s2"))
  expect_equal(names(x$phenotypes), "biomass")
})

test_that("mismatches are reported clearly", {
  o <- own()
  expect_error(pipeline_input(o$ph, c(5, 10), o$X), "Environment values not in environments")
  o$ph$biomass$Genotype[1] <- "var_Z"
  expect_error(pipeline_input(o$ph, c(5, 9), o$X), "genotypes not in snp_matrix: var_Z")
  o <- own()
  expect_error(pipeline_input(o$ph, c(5, 9), o$X, truth = data.frame(Shape = "x", Parameter = "p", SNP = "s1")),
               "truth\\$Shape must match")
  expect_error(pipeline_input(unname(o$ph), c(5, 9), o$X), "named list")
})

test_that("own data without SNPs is allowed; truth then is not", {
  o <- own()
  expect_null(pipeline_input(o$ph, c(5, 9))$snp_matrix)
  expect_error(pipeline_input(o$ph, c(5, 9), truth = data.frame(Shape = "biomass", Parameter = "p", SNP = "s1")),
               "truth needs a snp_matrix")
})
