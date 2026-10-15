# =============================================================================
# 06_analyze_gwas_results.R — summarise the GWAS results of all scenarios
# =============================================================================
# WHAT IT DOES: Summarises, per plasticity index, how many causal SNPs are found in
#   every heritability scenario (h2 x causal SNPs per parameter). A parameter is
#   always the parameter of one shape ("shape: parameter"); its causal SNPs are
#   searched only in the GWAS of that shape, and a hit on a SNP that is not causal
#   for that shape is a false positive of that shape. Draws count/fraction heatmaps,
#   precision and recall per PPI, true/false positives, and PPI coverage summaries.
# REQUIRES:     The tables of R/gwas/scenarios/combine_scenarios.R (run the scenarios first).
# PRODUCES:     Under COMBINED_DIR: plots/gwas_scenario_summary.pdf,
#               precision_recall_per_ppi.csv, recall_per_shape_parameter.csv.
# HOW TO RUN:   Rscript --vanilla R/gwas/06_analyze_gwas_results.R   (from the repo root)
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   COMBINED_DIR  path, default output/combined_scenarios   tables of combine_scenarios.R [COMMON]
#   PPIs          character, default all PPIs in the tables  which plasticity indices to summarise
#   FOCUS_PPIS    character, default c("RN", "RTR", "FW")    PPI set for the coverage summary
#   OTHER_PPIS    character, default c("APC", "D_slope", "PSI", "RNN")  second PPI set
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
library(pheatmap)
library(ggplot2)

if (!exists("COMBINED_DIR")) COMBINED_DIR <- here::here("output", "combined_scenarios")
if (!exists("FOCUS_PPIS")) FOCUS_PPIS <- c("RN", "RTR", "FW")
if (!exists("OTHER_PPIS")) OTHER_PPIS <- c("APC", "D_slope", "PSI", "RNN")

## ---- Import the combined tables ----
detections <- read.csv(file.path(COMBINED_DIR, "detection_counts.csv"), stringsAsFactors = FALSE)
hits       <- read.csv(file.path(COMBINED_DIR, "gwas_hits.csv"), stringsAsFactors = FALSE)
truth      <- read.csv(file.path(COMBINED_DIR, "causal_snps_truth.csv"), stringsAsFactors = FALSE)
p_threshold <- readRDS(file.path(COMBINED_DIR, "settings.rds"))$p_threshold

if (!exists("PPIs")) PPIs <- sort(unique(detections$Score))
detections <- detections[detections$Score %in% PPIs, ]
hits <- hits[hits$Score %in% PPIs, ]

## ---- Constants: shape-parameter pairs, scenarios ----
causal_rows <- detections[detections$Parameter != "False_Positives", ]
fp_rows     <- detections[detections$Parameter == "False_Positives", ]
RN_param_pairs <- unique(causal_rows[, c("Shape", "Parameter")])
RN_param_pairs <- RN_param_pairs[order(match(RN_param_pairs$Shape, c("linear", "gaussian", "sinusoidal", "wave")),
                                       RN_param_pairs$Parameter), ]
pair_labels <- paste0(RN_param_pairs$Shape, ": ", RN_param_pairs$Parameter)
shapes  <- unique(RN_param_pairs$Shape)
numSNPs <- sort(unique(detections$Mc))
H_vals  <- sort(unique(detections$H2))
H_pct   <- round(100 * H_vals)

## ---- Arrays: shape-parameter x PPI x Mc x h2 ----
dims <- c(nrow(RN_param_pairs), length(PPIs), length(numSNPs), length(H_vals))
foundCausalSNPsNum <- array(0, dim = dims)   # causal SNPs of this shape-parameter found
foundSNPsNum       <- array(0, dim = dims)   # all hits of this PPI in the GWAS of this shape
trueSNPsNum        <- array(0, dim = dims)   # causal SNPs of this shape-parameter
for (rnp_i in seq_len(nrow(RN_param_pairs))) {
  for (PPI_i in seq_along(PPIs)) {
    for (nSNP_i in seq_along(numSNPs)) {
      for (H_i in seq_along(H_vals)) {
        sel <- causal_rows$Shape == RN_param_pairs$Shape[rnp_i] & causal_rows$Parameter == RN_param_pairs$Parameter[rnp_i] &
          causal_rows$Score == PPIs[PPI_i] & causal_rows$Mc == numSNPs[nSNP_i] & causal_rows$H2 == H_vals[H_i]
        if (any(sel)) {
          foundCausalSNPsNum[rnp_i, PPI_i, nSNP_i, H_i] <- causal_rows$Found[sel]
          trueSNPsNum[rnp_i, PPI_i, nSNP_i, H_i] <- causal_rows$True_SNPs[sel]
        }
        hsel <- hits$Type == RN_param_pairs$Shape[rnp_i] & hits$Score == PPIs[PPI_i] &
          hits$Mc == numSNPs[nSNP_i] & hits$H2 == H_vals[H_i]
        foundSNPsNum[rnp_i, PPI_i, nSNP_i, H_i] <- length(unique(hits$SNP[hsel]))
      }
    }
  }
}

## ---- Sorted 2D matrices: rows shape-parameter x Mc, columns PPI x h2 ----
PPIs_sorted_i   <- order(apply(foundCausalSNPsNum, 2, sum), decreasing = TRUE)
RN_par_sorted_i <- order(apply(foundCausalSNPsNum, 1, sum), decreasing = TRUE)

to_2d <- function(A, row_order = seq_len(dims[1]), col_order = seq_len(dims[2])) {
  out <- matrix(0, dims[1] * dims[3], dims[2] * dims[4])
  for (r in seq_len(dims[1])) for (p in seq_len(dims[2])) for (n in seq_len(dims[3])) for (h in seq_len(dims[4])) {
    out[dims[3] * (r - 1) + n, dims[4] * (p - 1) + h] <- A[row_order[r], col_order[p], n, h]
  }
  out
}
row_labels <- function(order) as.vector(t(outer(pair_labels[order], paste(numSNPs, "SNPs"), paste)))
col_labels <- function(order) as.vector(t(outer(PPIs[order], paste0(H_pct, "% h2"), paste)))

fractions <- foundCausalSNPsNum / ifelse(trueSNPsNum > 0, trueSNPsNum, NA)

plots_dir <- file.path(COMBINED_DIR, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plots_dir, "gwas_scenario_summary.pdf"), width = 16, height = 12)

heat <- function(M, rl, cl, main) {
  pheatmap(M, labels_row = rl, labels_col = cl, cluster_rows = FALSE, cluster_cols = FALSE, main = main)
}
heat(to_2d(foundSNPsNum, RN_par_sorted_i, PPIs_sorted_i), row_labels(RN_par_sorted_i), col_labels(PPIs_sorted_i),
     paste("All SNPs with p <", p_threshold, "in the GWAS of the shape (sorted)"))
heat(to_2d(foundCausalSNPsNum, RN_par_sorted_i, PPIs_sorted_i), row_labels(RN_par_sorted_i), col_labels(PPIs_sorted_i),
     paste("Causal SNPs found, p <", p_threshold, "(sorted)"))
heat(to_2d(fractions, RN_par_sorted_i, PPIs_sorted_i), row_labels(RN_par_sorted_i), col_labels(PPIs_sorted_i),
     "Fraction of the causal SNPs found (sorted)")
heat(to_2d(foundSNPsNum), row_labels(seq_len(dims[1])), col_labels(seq_len(dims[2])),
     paste("All SNPs with p <", p_threshold, "in the GWAS of the shape"))
heat(to_2d(foundCausalSNPsNum), row_labels(seq_len(dims[1])), col_labels(seq_len(dims[2])),
     paste("Causal SNPs found, p <", p_threshold))

## ---- Histograms ----
hist(foundSNPsNum, main = paste("All SNPs with p <", p_threshold, "per scenario, shape-parameter and PPI"))
hist(foundCausalSNPsNum, main = paste("Causal SNPs with p <", p_threshold, "per scenario, shape-parameter and PPI"))

## ---- Causal SNPs found, summed over the scenarios ----
summedCausal <- rowSums(foundCausalSNPsNum, dims = 2)
bubble <- data.frame(PPI = factor(rep(PPIs, each = dims[1]), levels = PPIs),
                     Pair = factor(rep(pair_labels, times = dims[2]), levels = rev(pair_labels)),
                     Found = as.vector(summedCausal))
print(ggplot(bubble[bubble$Found > 0, ], aes(x = PPI, y = Pair, size = Found)) +
        geom_point(shape = 21, fill = "steelblue", alpha = 0.6) +
        scale_size_area(max_size = 10) +
        labs(x = NULL, y = NULL, title = "Causal SNPs found, summed over all scenarios") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))

## ---- True positives, false positives, precision and recall ----
# per shape-parameter x PPI, summed over the scenarios
truePositivesPerPPIScen <- rowSums(foundCausalSNPsNum, dims = 2)
allCausalsPerPPIScen    <- rowSums(trueSNPsNum, dims = 2)
# false positives are counted per shape (hits on SNPs that are not causal for the shape)
falsePositivesPerPPI_RN <- t(sapply(shapes, function(s) sapply(PPIs, function(p) {
  sum(fp_rows$Found[fp_rows$Shape == s & fp_rows$Score == p])
})))
if (length(PPIs) == 1) falsePositivesPerPPI_RN <- matrix(falsePositivesPerPPI_RN, ncol = 1)
dimnames(falsePositivesPerPPI_RN) <- list(shapes, PPIs)
truePositivesPerRN <- t(sapply(shapes, function(s) colSums(truePositivesPerPPIScen[RN_param_pairs$Shape == s, , drop = FALSE])))
if (length(PPIs) == 1) truePositivesPerRN <- matrix(truePositivesPerRN, ncol = 1)

truePositivesPerPPI  <- colSums(truePositivesPerPPIScen)
falsePositivesPerPPI <- colSums(falsePositivesPerPPI_RN)
allCausalsPerPPI     <- colSums(allCausalsPerPPIScen)
precisionPerPPI <- truePositivesPerPPI / (truePositivesPerPPI + falsePositivesPerPPI)
recallPerPPI    <- truePositivesPerPPI / allCausalsPerPPI
recallPerPPIScen <- truePositivesPerPPIScen / allCausalsPerPPIScen
dimnames(recallPerPPIScen) <- list(pair_labels, PPIs)
precisionPerPPI_RN <- truePositivesPerRN / (truePositivesPerRN + falsePositivesPerPPI_RN)
dimnames(precisionPerPPI_RN) <- list(shapes, PPIs)

write.csv(data.frame(PPI = PPIs, True_Positives = truePositivesPerPPI, False_Positives = falsePositivesPerPPI,
                     Causal_SNPs = allCausalsPerPPI, Precision = precisionPerPPI, Recall = recallPerPPI),
          file.path(COMBINED_DIR, "precision_recall_per_ppi.csv"), row.names = FALSE)
write.csv(data.frame(Shape = rep(RN_param_pairs$Shape, times = dims[2]), Parameter = rep(RN_param_pairs$Parameter, times = dims[2]),
                     PPI = rep(PPIs, each = dims[1]), Found = as.vector(truePositivesPerPPIScen),
                     Causal_SNPs = as.vector(allCausalsPerPPIScen), Recall = as.vector(recallPerPPIScen)),
          file.path(COMBINED_DIR, "recall_per_shape_parameter.csv"), row.names = FALSE)

# precision (left axis) and recall (right axis, scaled) per PPI
sf <- max(precisionPerPPI, na.rm = TRUE) / max(max(recallPerPPI, na.rm = TRUE), 1e-12)
val_df <- data.frame(PPI = rep(PPIs, 2), Value = c(precisionPerPPI, recallPerPPI * sf),
                     Measure = rep(c("precision", "recall"), each = length(PPIs)))
print(ggplot(val_df, aes(x = PPI, y = Value, fill = Measure)) +
        geom_col(position = "dodge") +
        scale_y_continuous(sec.axis = sec_axis(~ . / sf, name = "recall")) +
        labs(y = "precision", title = "Precision and recall per PPI (all shapes and scenarios)") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))

# recall per PPI, one panel per shape-parameter
rec_df <- data.frame(PPI = rep(PPIs, each = dims[1]), Pair = factor(rep(pair_labels, times = dims[2]), levels = pair_labels),
                     Recall = as.vector(recallPerPPIScen))
print(ggplot(rec_df, aes(x = PPI, y = Recall)) + geom_col(fill = "steelblue") +
        facet_wrap(~Pair) + labs(title = "Recall per PPI, per shape-parameter") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)))

pheatmap(recallPerPPIScen, cluster_rows = FALSE, cluster_cols = FALSE, main = "Recall per shape-parameter and PPI")
pheatmap(precisionPerPPI_RN, cluster_rows = FALSE, cluster_cols = FALSE, main = "Precision per shape and PPI")
dimnames(truePositivesPerPPIScen) <- list(pair_labels, PPIs)
pheatmap(truePositivesPerPPIScen, cluster_rows = FALSE, cluster_cols = FALSE, main = "True positives per shape-parameter and PPI")
pheatmap(falsePositivesPerPPI_RN, cluster_rows = FALSE, cluster_cols = FALSE, main = "False positives per shape and PPI")

# true positives of every shape-parameter, followed by the false positives of that shape (log scale)
tp_fp <- do.call(rbind, lapply(shapes, function(s) {
  rbind(truePositivesPerPPIScen[RN_param_pairs$Shape == s, , drop = FALSE],
        matrix(falsePositivesPerPPI_RN[s, ], nrow = 1, dimnames = list(paste(s, "false positives"), PPIs)))
}))
log_tp_fp <- log(tp_fp)
log_tp_fp[!is.finite(log_tp_fp) | log_tp_fp < 0] <- NA
breaks <- c(1, 4, 10, 25, 50, 100, 200)
if (any(is.finite(log_tp_fp))) {
  pheatmap(log_tp_fp, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 12,
           gaps_row = cumsum(table(factor(RN_param_pairs$Shape, levels = shapes)) + 1)[-length(shapes)],
           angle_col = 90, legend_breaks = log(breaks), legend_labels = breaks,
           main = "True positives per shape-parameter and false positives per shape (log scale)")
}
dev.off()

## ---- PPI coverage: which causal SNPs are found by which PPIs ----
# one row per found causal SNP (scenario x shape-parameter x SNP), one column per PPI
causal_hits <- merge(hits, truth, by.x = c("H2", "Mc", "Type", "SNP"), by.y = c("H2", "Mc", "Shape", "SNP"))
causal_hits$Key <- paste(causal_hits$H2, causal_hits$Mc, causal_hits$Type, causal_hits$Parameter, causal_hits$SNP)
keys <- unique(causal_hits[, c("Key", "Type", "Parameter")])
ppiFoundCausalSNPsMatrix <- matrix(FALSE, nrow(keys), length(PPIs), dimnames = list(keys$Key, PPIs))
ppiFoundCausalSNPsMatrix[cbind(match(causal_hits$Key, keys$Key), match(causal_hits$Score, PPIs))] <- TRUE

found_by <- function(set, rows = rep(TRUE, nrow(keys))) {
  set <- intersect(set, PPIs)
  sum(rowSums(ppiFoundCausalSNPsMatrix[rows, set, drop = FALSE]) > 0)
}
total_found_causals <- nrow(keys)
message("Number of all found causal SNPs: ", total_found_causals)
message("  found by ", paste(FOCUS_PPIS, collapse = ", "), ": ", found_by(FOCUS_PPIS))
message("  found by ", paste(OTHER_PPIS, collapse = ", "), ": ", found_by(OTHER_PPIS))
for (i in seq_len(nrow(RN_param_pairs))) {
  rows <- keys$Type == RN_param_pairs$Shape[i] & keys$Parameter == RN_param_pairs$Parameter[i]
  if (!any(rows)) next
  message("  ", pair_labels[i], ": ", sum(rows), " found; by ", paste(FOCUS_PPIS, collapse = ", "), ": ",
          found_by(FOCUS_PPIS, rows), " (", round(100 * found_by(FOCUS_PPIS, rows) / sum(rows)), " %)")
}
message("Causal SNPs found by none of ", paste(FOCUS_PPIS, collapse = ", "), ", per PPI:")
print(colSums(ppiFoundCausalSNPsMatrix[rowSums(ppiFoundCausalSNPsMatrix[, intersect(FOCUS_PPIS, PPIs), drop = FALSE]) == 0, , drop = FALSE]))

message("Plots  : ", file.path(plots_dir, "gwas_scenario_summary.pdf"))
message("Tables : ", file.path(COMBINED_DIR, "precision_recall_per_ppi.csv"), ", recall_per_shape_parameter.csv")
