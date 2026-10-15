# =============================================================================
# 06_analyze_maize_gwas_results.R — summary of the maize GWAS
# =============================================================================
# WHAT IT DOES: Summarises the GWAS hits of the maize data per plasticity index (PPI)
#   and trait (biomass, leaf area, WUE): which SNPs each PPI finds, how many hits per
#   PPI and trait, the allele frequency of the hit SNPs, and how similar the PPIs are.
#   Hits: significant = Bonferroni P < ALPHA / number of SNPs; additionally all SNPs
#   with -log10(P) > 3 and > 5.
# REQUIRES:     The maize resolution folder with PPIs and GWAS (scenario_maize.R).
# PRODUCES:     In RESOLUTION_DIR/<REPLICATE_HANDLING>/gwas/: maize_hits_significant.csv, maize_hits_logp3.csv,
#               maize_hit_counts.csv, maize_hit_snps.csv (SNP x PPI per trait, allele
#               frequency); RESOLUTION_DIR/<REPLICATE_HANDLING>/plots/maize_gwas_summary.pdf
# HOW TO RUN:   setwd("~/PlastQuant")
#               RESOLUTION_DIR <- here::here("output", "populations", "maize", "n_8")
#               source(here::here("R", "gwas", "06_analyze_maize_gwas_results.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   RESOLUTION_DIR  path, default output/populations/maize/n_8   maize resolution folder [COMMON]
#   REPLICATE_HANDLING  "mean_before" (default) or "replicates"   which PPIs of the resolution
#   ALPHA           numeric, default 0.05                         Bonferroni level
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded
library(pheatmap)
source(here::here("R", "pipeline", "load_pipeline.R"))

if (!exists("RESOLUTION_DIR")) RESOLUTION_DIR <- here::here("output", "populations", "maize", "n_8")
if (!exists("REPLICATE_HANDLING")) REPLICATE_HANDLING <- "mean_before"   # PPIs of mean_before/ or replicates/
PPI_SET_DIR <- file.path(RESOLUTION_DIR, REPLICATE_HANDLING)
if (!exists("ALPHA")) ALPHA <- 0.05

## ---- GWAS results and PPIs of the maize data ----
gwas <- read.csv(file.path(PPI_SET_DIR, "gwas", "gwas_results.csv"), stringsAsFactors = FALSE)
n_snps <- length(unique(gwas$SNP))
threshold <- ALPHA / n_snps
significant <- gwas[gwas$P_Value < threshold, ]
logp3 <- gwas[gwas$LogP > 3, ]
logp5 <- gwas[gwas$LogP > 5, ]
message("Maize GWAS: ", n_snps, " SNPs; ", nrow(significant), " significant hits (P < ",
        format(threshold, digits = 2), "), ", nrow(logp3), " with -log10 P > 3, ", nrow(logp5), " > 5")

PPIs <- sort(unique(gwas$Score))
RNs  <- c("leafArea", "WUE", "biomass")

## ---- Hit counts per PPI and trait ----
count_table <- function(hits, label) {
  tab <- table(factor(hits$Score, levels = PPIs), factor(hits$Type, levels = RNs))
  data.frame(Hits = label, PPI = rep(PPIs, times = length(RNs)), Trait = rep(RNs, each = length(PPIs)),
             N_SNPs = as.vector(tab), stringsAsFactors = FALSE)
}
counts <- rbind(count_table(significant, "significant"), count_table(logp3, "logp3"), count_table(logp5, "logp5"))

## ---- Which PPI finds which SNP (significant hits), per trait and over all traits ----
SNPs <- sort(unique(significant$SNP))
ppi_found <- function(hits) {
  m <- matrix(FALSE, length(SNPs), length(PPIs), dimnames = list(SNPs, PPIs))
  if (nrow(hits)) m[cbind(match(hits$SNP, SNPs), match(hits$Score, PPIs))] <- TRUE
  m
}
ppiFoundSNPsMatrix <- ppi_found(significant)
ppiFoundSNPsMatrixPerTrait <- lapply(setNames(RNs, RNs), function(rn) ppi_found(significant[significant$Type == rn, ]))

## ---- Allele frequency of the hit SNPs ----
pop <- read_population(RESOLUTION_DIR)
X <- pop$snp_matrix[, SNPs, drop = FALSE]
p <- colMeans(X) / 2
hit_snps <- data.frame(SNP = SNPs, MAF = pmin(p, 1 - p), Carrier_Frequency = colMeans(X > 0),
                       N_PPIs = rowSums(ppiFoundSNPsMatrix), stringsAsFactors = FALSE)
for (rn in RNs) hit_snps[[paste0("PPIs_", rn)]] <- rowSums(ppiFoundSNPsMatrixPerTrait[[rn]])

## ---- Similarity of the PPIs (pooled over the three traits) ----
traits <- readRDS(file.path(PPI_SET_DIR, "ppis", "ppi_scores.rds"))
scores <- ppi_wide_table(traits, shapes = RNs)
cor_scores <- cor(scores[, setdiff(names(scores), c("Shape", "Genotype"))], use = "pairwise.complete.obs",
                  method = "kendall")

## ---- Save ----
gwas_dir <- population_subdir(PPI_SET_DIR, "gwas")
write.csv(significant, file.path(gwas_dir, "maize_hits_significant.csv"), row.names = FALSE)
write.csv(logp3, file.path(gwas_dir, "maize_hits_logp3.csv"), row.names = FALSE)
write.csv(counts, file.path(gwas_dir, "maize_hit_counts.csv"), row.names = FALSE)
write.csv(hit_snps, file.path(gwas_dir, "maize_hit_snps.csv"), row.names = FALSE)

plots_dir <- population_subdir(PPI_SET_DIR, "plots")
pdf(file.path(plots_dir, "maize_gwas_summary.pdf"), width = 12, height = 8)
for (h in c("significant", "logp3", "logp5")) {
  m <- matrix(counts$N_SNPs[counts$Hits == h], length(PPIs), length(RNs), dimnames = list(PPIs, RNs))
  pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%d",
           main = paste0("Maize: SNP hits per PPI and trait (", h, ")"))
}
if (length(SNPs) > 1) {
  pheatmap(ppiFoundSNPsMatrix * 1, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
           main = "Significant SNPs (rows) found by each PPI (columns), all traits")
}
pheatmap(cor_scores, main = "Kendall correlation of the PPIs (maize, all traits)")
dev.off()
message("Maize GWAS summary written to ", gwas_dir, " and ", plots_dir)
