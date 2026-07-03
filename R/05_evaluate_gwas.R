# Script: evaluate_gwas_accuracy.R
# Purpose: Compare GWAS results against simulated ground truth to validate plasticity scores.
#          Produces 5 power heatmaps: one combined (all forms) + one per RN form.
source(here::here("R", "04_run_gwas.R"))
library(dplyr)
library(ggplot2)
library(tidyr)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
results_dir <- latest_run
plots_dir   <- file.path(latest_run, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

gwas_file <- file.path(results_dir, "gwas_results_mlm_80.csv")
if (!file.exists(gwas_file)) gwas_file <- file.path(results_dir, "gwas_results_lm_80.csv")
if (!file.exists(gwas_file)) stop("GWAS results file not found. Please run run_gwas_simulation.R first.")

message("Loading GWAS results from: ", gwas_file)
gwas_res <- read.csv(gwas_file)

truth_file <- file.path(results_dir, "causal_snps_truth.csv")
message("Loading truth from: ", truth_file)
truth <- read.csv(truth_file)

# -----------------------------------------------------------------------------
# 2. Significance Threshold & Metadata
# -----------------------------------------------------------------------------
alpha   <- 0.05
n_snps  <- length(unique(gwas_res$SNP))
threshold <- alpha / n_snps
message("Significance threshold (Bonferroni): P < ", format(threshold, scientific = TRUE))

if (exists("snp_matrix") && is.matrix(snp_matrix)) {
  n_genotypes <- nrow(snp_matrix)
} else if (exists("kinship_matrix") && is.matrix(kinship_matrix)) {
  n_genotypes <- nrow(kinship_matrix)
} else {
  n_genotypes <- NA
}

subtitle_base <- paste0(
  "N_Genotypes = ", n_genotypes,
  ", N_SNPs = ", n_snps,
  "  (Bonferroni P < ", format(threshold, digits = 2, scientific = TRUE), ")"
)

# -----------------------------------------------------------------------------
# 3. Helper: compute power table for a given set of significant hits
# -----------------------------------------------------------------------------
compute_power <- function(sig_hits_subset, truth, scores, params) {
  rows <- list()
  for (s in scores) {
    score_hits <- unique(sig_hits_subset$SNP[sig_hits_subset$Score == s])
    for (p in params) {
      true_snps <- truth$SNP[truth$Parameter == p]
      n_true    <- length(true_snps)
      n_found   <- length(intersect(score_hits, true_snps))
      rows[[paste(s, p)]] <- data.frame(
        Score     = s,
        Parameter = p,
        True_SNPs = n_true,
        Found     = n_found,
        Power     = if (n_true > 0) n_found / n_true else NA_real_
      )
    }
    # Global false positives
    n_fp <- length(setdiff(score_hits, unique(truth$SNP)))
    rows[[paste(s, "FP")]] <- data.frame(
      Score     = s,
      Parameter = "False_Positives_Global",
      True_SNPs = NA_integer_,
      Found     = n_fp,
      Power     = NA_real_
    )
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# 4. Helper: draw power heatmap
# -----------------------------------------------------------------------------
make_heatmap <- function(power_df, title, subtitle) {
  df <- power_df %>%
    filter(Parameter != "False_Positives_Global") %>%
    mutate(Power = as.numeric(Power))

  ggplot(df, aes(x = Parameter, y = Score, fill = Power)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), na.value = "grey90") +
    geom_text(aes(label = ifelse(is.na(Power), "", sprintf("%.2f", Power))),
              color = "black", size = 3) +
    labs(title = title, subtitle = subtitle,
         x = "Causal Parameter (Truth)", y = "Plasticity Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------------------------------------------------------
# 5. Evaluate: combined (union of hits across all forms) + per form
# -----------------------------------------------------------------------------
sig_hits <- gwas_res %>% filter(P_Value < threshold)
scores   <- unique(gwas_res$Score)
params   <- unique(truth$Parameter)
forms    <- unique(gwas_res$Type)

# Combined evaluation (union of SNP hits across all forms per score)
eval_combined <- compute_power(sig_hits, truth, scores, params)
eval_combined$Form <- "All forms (union)"

# Per-form evaluation
eval_per_form <- lapply(forms, function(f) {
  hits_f <- sig_hits %>% filter(Type == f)
  df     <- compute_power(hits_f, truth, scores, params)
  df$Form <- f
  df
})
names(eval_per_form) <- forms

# -----------------------------------------------------------------------------
# 6. Save CSVs
# -----------------------------------------------------------------------------
all_eval <- bind_rows(eval_combined, bind_rows(eval_per_form))
write.csv(all_eval,     file.path(results_dir, "gwas_accuracy_summary.csv"),          row.names = FALSE)
write.csv(eval_combined, file.path(results_dir, "gwas_accuracy_summary_combined.csv"), row.names = FALSE)
for (f in forms) {
  write.csv(eval_per_form[[f]],
            file.path(results_dir, paste0("gwas_accuracy_summary_", f, ".csv")),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 7. Plot: 5 heatmaps in one PDF
# -----------------------------------------------------------------------------
pdf(file.path(plots_dir, "gwas_power_heatmap.pdf"), width = 10, height = 8)

# Page 1: combined
print(make_heatmap(eval_combined,
                   "GWAS Power — All Forms (union of hits)",
                   subtitle_base))

# Pages 2–5: one per form
for (f in forms) {
  n_form <- nrow(sig_hits %>% filter(Type == f) %>% distinct(SNP, Score))
  print(make_heatmap(eval_per_form[[f]],
                     paste0("GWAS Power — Form: ", f),
                     paste0(subtitle_base, "  |  GWAS run on ", f, " genotypes only")))
}

dev.off()

message("Evaluation complete.")
message("Summary CSV : ", file.path(results_dir, "gwas_accuracy_summary.csv"))
message("Heatmap PDF : ", file.path(plots_dir, "gwas_power_heatmap.pdf"))
