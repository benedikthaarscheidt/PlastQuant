# =============================================================================
# pipeline_input.R — the data every pipeline run works on, simulated or own
# =============================================================================
# A pipeline input is a list with
#   snp_matrix    genotypes x SNPs, dosages 0/1/2, rownames = genotype IDs;
#                 optional for own data (then only the PPIs are computed, no GWAS)
#   environments  numeric vector of environment values (one per measurement point)
#   phenotypes    named list, one data.frame per shape or trait with columns
#                 Genotype, Replicate, Environment, Trait
#   truth         optional data.frame Shape, Parameter, SNP (known causal SNPs);
#                 Shape must match a name of `phenotypes`
# simulate_reaction_norms() returns exactly this (plus simulation details).
# =============================================================================

pipeline_input <- function(phenotypes, environments, snp_matrix = NULL, truth = NULL) {
  x <- list(snp_matrix = snp_matrix, phenotypes = phenotypes, environments = environments, truth = truth)
  check_pipeline_input(x)
  x
}

check_pipeline_input <- function(x) {
  has_snps <- !is.null(x$snp_matrix)
  if (has_snps) {
    if (!is.matrix(x$snp_matrix) || is.null(rownames(x$snp_matrix)) || is.null(colnames(x$snp_matrix))) {
      stop("snp_matrix must be a matrix with genotype IDs as rownames and SNP IDs as colnames")
    }
    if (!all(x$snp_matrix %in% c(0, 1, 2))) stop("snp_matrix must contain dosages 0, 1, 2")
  }
  if (!is.numeric(x$environments) || length(x$environments) < 2) stop("environments must be numeric, length >= 2")
  if (!is.list(x$phenotypes) || is.null(names(x$phenotypes)) || any(names(x$phenotypes) == "")) {
    stop("phenotypes must be a named list (one table per shape or trait)")
  }
  for (s in names(x$phenotypes)) {
    d <- x$phenotypes[[s]]
    miss <- setdiff(c("Genotype", "Replicate", "Environment", "Trait"), names(d))
    if (length(miss)) stop("phenotypes$", s, " lacks column(s): ", paste(miss, collapse = ", "))
    unknown <- if (has_snps) setdiff(unique(as.character(d$Genotype)), rownames(x$snp_matrix)) else character(0)
    if (length(unknown)) stop("phenotypes$", s, ": genotypes not in snp_matrix: ", paste(head(unknown, 5), collapse = ", "))
    bad_env <- setdiff(unique(d$Environment), x$environments)
    if (length(bad_env)) stop("phenotypes$", s, ": Environment values not in environments: ", paste(head(bad_env, 5), collapse = ", "))
  }
  if (!is.null(x$truth)) {
    if (!has_snps) stop("truth needs a snp_matrix")
    miss <- setdiff(c("Shape", "Parameter", "SNP"), names(x$truth))
    if (length(miss)) stop("truth lacks column(s): ", paste(miss, collapse = ", "))
    if (!all(x$truth$SNP %in% colnames(x$snp_matrix))) stop("truth contains SNPs that are not in snp_matrix")
    if (!all(x$truth$Shape %in% names(x$phenotypes))) stop("truth$Shape must match the names of phenotypes")
  }
  invisible(TRUE)
}
