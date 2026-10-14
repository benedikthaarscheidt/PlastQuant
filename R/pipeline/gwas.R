# =============================================================================
# gwas.R — GWAS of every PPI, separately for every shape (rrBLUP mixed model)
# =============================================================================

# Realised relationship matrix from 0/1/2 dosages
kinship_matrix <- function(snp_matrix) rrBLUP::A.mat(snp_matrix - 1)

# One GWAS: trait values y named by line (genotype) ID against all SNPs. A line may
# occur several times (one value per replicate); all its rows get the SNPs of the line,
# and the kinship links them, so replicates of a line are not counted as separate lines.
gwas_one <- function(y, snp_matrix, kinship) {
  y <- y[is.finite(y)]
  ids <- names(y)
  lines <- unique(ids)
  geno <- data.frame(name = colnames(snp_matrix), chrom = 1, pos = seq_len(ncol(snp_matrix)),
                     t(snp_matrix[lines, , drop = FALSE] - 1), check.names = FALSE)
  # rrBLUP prints progress lines; they are captured to keep the console readable
  utils::capture.output(
    res <- rrBLUP::GWAS(data.frame(line = ids, y = unname(y)), geno, K = kinship[lines, lines],
                        min.MAF = 0, n.PC = 0, P3D = TRUE, plot = FALSE)
  )
  data.frame(SNP = res$name, LogP = res$y, P_Value = 10^(-res$y), stringsAsFactors = FALSE)
}

# traits: the PPI values used as GWAS traits (output of compute_ppis());
# returns one table with the columns SNP, LogP, P_Value, Score (PPI) and Type (shape)
run_gwas <- function(traits, snp_matrix, kinship = kinship_matrix(snp_matrix)) {
  out <- list()
  for (score_name in names(traits)) {
    sc <- traits[[score_name]]
    for (shape in unique(sc$Type)) {
      d <- sc[sc$Type == shape, ]
      y <- setNames(d$Score, d$Genotype_ID)
      if (sum(is.finite(y)) < 3) {
        message("GWAS skipped for ", score_name, " / ", shape, ": fewer than 3 values")
        next
      }
      r <- gwas_one(y, snp_matrix, kinship)
      r$Score <- score_name
      r$Type <- shape
      out[[length(out) + 1]] <- r
    }
  }
  do.call(rbind, out)
}
