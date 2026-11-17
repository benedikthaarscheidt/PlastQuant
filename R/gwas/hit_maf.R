# =============================================================================
# hit_maf.R — true parameters as GWAS traits; minor allele frequency of the GWAS hits
# =============================================================================
# Parts of run_population_gwas():
#   true_parameter_traits()  the parameters of the genotypic (mother) reaction norms as
#                            traits next to the PPIs, the "perfect" traits for the GWAS
#   hit_maf_tables()         MAF of the significant hits (true positive, in LD with a causal
#                            SNP, false positive) and of the causal SNPs (found by how many traits)
#   plot_hit_maf()           both as figures

# Minor allele frequency of every SNP (column) of a 0/1/2 dosage matrix
snp_maf <- function(X) {
  f <- colMeans(X) / 2
  pmin(f, 1 - f)
}

# The parameters of the mother reaction norms (data/mother_parameters.csv: Genotype, Shape,
# one column per parameter) as traits in the format of the PPIs (Score, Type, Genotype_ID,
# Replicate), one per parameter, named true_<parameter>. The stability V enters through
# the repeatability rho of the genotype, a decreasing function of its genetic value for V
# (genetic_values.R): true_V. Parameters a shape does not have are left out for that shape.
true_parameter_traits <- function(mother_parameters) {
  mp <- mother_parameters
  cols <- setdiff(names(mp), c("Genotype", "Shape"))
  out <- lapply(cols, function(p) {
    d <- mp[is.finite(mp[[p]]), , drop = FALSE]
    data.frame(Score = d[[p]], Type = d$Shape, Genotype_ID = as.character(d$Genotype),
               Replicate = rep(NA_integer_, nrow(d)), stringsAsFactors = FALSE)
  })
  names(out) <- paste0("true_", ifelse(cols == "rho", "V", cols))
  out[vapply(out, nrow, integer(1)) > 0]
}

# For every causal SNP of one shape the SNPs whose hit finds it, with the rules of
# compute_power(): its tag SNPs on a SNP panel (tags: SNP, Tag), the SNPs within `window`
# of it with LD (snp_pos: position of every SNP), otherwise the causal SNP itself
causal_hit_sets <- function(truth_shape, tags = NULL, window = 0, snp_pos = NULL) {
  cs <- unique(truth_shape$SNP)
  sets <- if (!is.null(tags)) {
    ts <- split(tags$Tag, tags$SNP)
    lapply(cs, function(s) unique(c(s, ts[[s]])))
  } else if (window > 0) {
    lapply(cs, function(s) names(snp_pos)[abs(snp_pos - snp_pos[s]) <= window])
  } else {
    as.list(cs)
  }
  setNames(sets, cs)
}

# gwas_results: run_gwas() output; maf: snp_maf() of the tested SNPs; threshold: P value
# below which a SNP is a hit. Returns
#   hits    one row per hit (Score, Type, SNP, P_Value, MAF, Status); with a truth table
#           Status is "true positive" (a causal SNP of its shape), "in LD with a causal
#           SNP" (one of its tag SNPs, or within the LD window) or "false positive",
#           without one "positive"
#   causal  (with truth) one row per causal SNP and shape: Parameter, MAF and the number
#           of PPIs and of true parameters (true_*) that found it
hit_maf_tables <- function(gwas_results, maf, threshold, truth = NULL, tags = NULL, window = 0, snp_pos = NULL) {
  sig <- gwas_results[gwas_results$P_Value < threshold, c("Score", "Type", "SNP", "P_Value"), drop = FALSE]
  sig$MAF <- unname(maf[sig$SNP])
  sig$Status <- rep(if (is.null(truth)) "positive" else "false positive", nrow(sig))
  causal <- NULL
  if (!is.null(truth)) {
    rows <- list()
    for (shape in unique(gwas_results$Type)) {
      tr <- truth[truth$Shape == shape, , drop = FALSE]
      if (!nrow(tr)) next
      sets <- causal_hit_sets(tr, tags, window, snp_pos)
      in_shape <- sig$Type == shape
      sig$Status[in_shape & sig$SNP %in% unlist(sets)] <- "in LD with a causal SNP"
      sig$Status[in_shape & sig$SNP %in% tr$SNP] <- "true positive"
      for (k in seq_len(nrow(tr))) {
        by <- unique(sig$Score[in_shape & sig$SNP %in% sets[[tr$SNP[k]]]])
        rows[[length(rows) + 1]] <- data.frame(
          Shape = shape, Parameter = tr$Parameter[k], SNP = tr$SNP[k], MAF = unname(maf[tr$SNP[k]]),
          Found_by_PPIs = sum(!grepl("^true_", by)), Found_by_true_parameters = sum(grepl("^true_", by)),
          stringsAsFactors = FALSE)
      }
    }
    causal <- do.call(rbind, rows)
  }
  list(hits = sig, causal = causal)
}

# Page 1, per shape, hits of the PPIs only (not of the true parameters), every SNP once per
# shape: MAF classes of all tested SNPs (reference), the positives (SNPs hit by any PPI),
# split with a truth table into true positives (a causal SNP), SNPs in LD with a causal SNP
# (neither true nor false) and false positives, and the false negatives (causal SNPs that
# no PPI found, neither on the SNP nor on a SNP in LD with it). Are the SNPs found rarer or
# more common than the SNPs tested?
# Page 2 (with truth): MAF of every causal SNP against the number of PPIs that found it;
# filled points: also found by a true parameter.
plot_hit_maf <- function(tabs, maf, file, subtitle = "") {
  breaks <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  maf_class <- function(x) cut(x, breaks, include.lowest = TRUE)
  share <- function(x) { t <- table(maf_class(x)); data.frame(MAF_class = names(t), Share = as.vector(t) / max(sum(t), 1)) }
  ppi_hits <- tabs$hits[!grepl("^true_", tabs$hits$Score), , drop = FALSE]
  hits <- unique(ppi_hits[, c("Type", "SNP", "MAF", "Status")])
  with_truth <- !is.null(tabs$causal)
  shapes <- unique(c(tabs$hits$Type, if (with_truth) tabs$causal$Shape))
  if (!length(shapes)) shapes <- "no hits"
  bars <- do.call(rbind, lapply(shapes, function(s) {
    h <- hits[hits$Type == s, , drop = FALSE]
    groups <- list(all = maf, positive = h$MAF)
    panel <- sprintf("%s: %d positive", s, nrow(h))
    if (with_truth) {
      fn <- tabs$causal$MAF[tabs$causal$Shape == s & tabs$causal$Found_by_PPIs == 0]
      groups <- c(groups, list(`true positive` = h$MAF[h$Status == "true positive"],
                               `in LD with a causal SNP` = h$MAF[h$Status == "in LD with a causal SNP"],
                               `false positive` = h$MAF[h$Status == "false positive"], `false negative` = fn))
      panel <- sprintf("%s (%d true, %d in LD, %d false), %d false negative", panel,
                       length(groups$`true positive`), length(groups$`in LD with a causal SNP`),
                       length(groups$`false positive`), length(fn))
    }
    do.call(rbind, lapply(names(groups), function(g) cbind(share(groups[[g]]), Kind = g, Shape = panel)))
  }))
  used <- levels(maf_class(maf))[levels(maf_class(maf)) %in% bars$MAF_class[bars$Share > 0]]
  bars <- bars[bars$MAF_class %in% used, ]
  bars$MAF_class <- factor(bars$MAF_class, levels = used)
  cols <- c("all" = "grey75", "positive" = "grey30", "true positive" = "#2C7BB6",
            "in LD with a causal SNP" = "#ABD9E9", "false positive" = "#D7191C", "false negative" = "#FDAE61")
  bars$Kind <- factor(bars$Kind, levels = names(cols))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 11, height = 7)
  on.exit(grDevices::dev.off())
  print(ggplot2::ggplot(bars, ggplot2::aes(x = MAF_class, y = Share, fill = Kind)) +
          ggplot2::geom_col(position = ggplot2::position_dodge(preserve = "single"), width = 0.85) +
          ggplot2::scale_fill_manual(values = cols, breaks = names(cols), name = NULL,
                                     labels = c("all tested SNPs", "positive", "true positive", "in LD with a causal SNP",
                                                "false positive", "false negative")) +
          ggplot2::facet_wrap(~ Shape) +
          ggplot2::labs(x = "Minor allele frequency", y = "Share of the SNPs of the group",
                        title = "MAF of the SNPs found by the PPIs, against all tested SNPs",
                        subtitle = paste0(subtitle, "\nPPIs only; every SNP counted once per shape. True positive: a causal SNP; in LD: r^2 >= 0.1 with a causal SNP;\n",
                                          "false negative: a causal SNP that no PPI found, neither itself nor in LD")) +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
  if (!is.null(tabs$causal) && nrow(tabs$causal)) {
    cz <- tabs$causal
    cz$Found_by_true <- cz$Found_by_true_parameters > 0
    print(ggplot2::ggplot(cz, ggplot2::aes(x = MAF, y = Found_by_PPIs, colour = Parameter, shape = Found_by_true)) +
            ggplot2::geom_point(position = ggplot2::position_jitter(width = 0, height = 0.15, seed = 1), size = 2.2, stroke = 0.7) +
            ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "Found by a true parameter") +
            ggplot2::facet_wrap(~ Shape) +
            ggplot2::labs(x = "Minor allele frequency of the causal SNP", y = "Number of PPIs that found it",
                          title = "Causal SNPs: minor allele frequency and detection", subtitle = subtitle) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::theme(legend.position = "top"))
  }
  invisible(file)
}
