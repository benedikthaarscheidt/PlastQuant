# =============================================================================
# gwas.R — GWAS of every PPI, separately for every shape (rrBLUP mixed model)
# =============================================================================

# Realised relationship matrix from 0/1/2 dosages
kinship_matrix <- function(snp_matrix) rrBLUP::A.mat(snp_matrix - 1)

# GWAS of several traits with the same rows in one rrBLUP call, so the SNP table is
# built only once. pheno: line (genotype) ID in column 1, one column per trait, NA where
# a trait has no value. A line may occur several times (one row per replicate); all its
# rows get the SNPs of the line, and the kinship links them, so replicates of a line are
# not counted as separate lines. rrBLUP fits every trait on its own rows with a value
# (lines and kinship subset to them), as a call with that trait alone would.
# n_cores > 1 splits the SNP tests over that many forked processes (same results).
gwas_traits <- function(pheno, snp_matrix, kinship, n_cores = 1) {
  lines <- unique(pheno[[1]])
  geno <- data.frame(name = colnames(snp_matrix), chrom = 1, pos = seq_len(ncol(snp_matrix)),
                     t(snp_matrix[lines, , drop = FALSE] - 1), check.names = FALSE)
  # rrBLUP prints progress lines; they are captured to keep the console readable
  utils::capture.output(
    res <- rrBLUP::GWAS(pheno, geno, K = kinship[lines, lines],
                        min.MAF = 0, n.PC = 0, n.core = n_cores, P3D = TRUE, plot = FALSE)
  )
  lapply(setNames(names(pheno)[-1], names(pheno)[-1]), function(tr) {
    data.frame(SNP = res$name, LogP = res[[tr]], P_Value = 10^(-res[[tr]]), stringsAsFactors = FALSE)
  })
}

# traits: the PPI values used as GWAS traits (output of compute_ppis());
# returns one table with the columns SNP, LogP, P_Value, Score (PPI) and Type (shape).
# Per shape, the PPIs with the same rows (genotype and replicate) share one rrBLUP call.
# n_cores: processes for the SNP tests (see gwas_traits())
run_gwas <- function(traits, snp_matrix, kinship = kinship_matrix(snp_matrix), n_cores = 1) {
  res <- list()   # res[[score]][[shape]]
  shapes <- unique(unlist(lapply(traits, function(sc) unique(sc$Type))))
  for (shape in shapes) {
    groups <- list()   # PPIs with identical rows -> one rrBLUP call
    for (score_name in names(traits)) {
      sc <- traits[[score_name]]
      d <- sc[sc$Type == shape, ]
      if (nrow(d) == 0) next
      if (sum(is.finite(d$Score)) < 3) {
        message("GWAS skipped for ", score_name, " / ", shape, ": fewer than 3 values")
        next
      }
      rep_id <- if (is.null(d$Replicate)) rep(NA, nrow(d)) else d$Replicate
      key <- paste(paste(d$Genotype_ID, rep_id, sep = "\r"), collapse = "\n")
      if (is.null(groups[[key]])) groups[[key]] <- list(line = d$Genotype_ID, scores = list())
      groups[[key]]$scores[[score_name]] <- ifelse(is.finite(d$Score), d$Score, NA_real_)
    }
    for (g in groups) {
      ids <- paste0("trait_", seq_along(g$scores))   # syntactic column names for rrBLUP
      pheno <- data.frame(line = g$line, setNames(g$scores, ids), stringsAsFactors = FALSE)
      r <- gwas_traits(pheno, snp_matrix, kinship, n_cores)
      for (k in seq_along(ids)) res[[names(g$scores)[k]]][[shape]] <- r[[ids[k]]]
    }
  }
  # rows in the order PPI, then shape, as before
  out <- list()
  for (score_name in names(traits)) {
    for (shape in unique(traits[[score_name]]$Type)) {
      r <- res[[score_name]][[shape]]
      if (is.null(r)) next
      r$Score <- score_name
      r$Type <- shape
      out[[length(out) + 1]] <- r
    }
  }
  do.call(rbind, out)
}
