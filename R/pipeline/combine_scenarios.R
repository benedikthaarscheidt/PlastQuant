# =============================================================================
# combine_scenarios.R — collect the populations of several scenarios into one table set
# =============================================================================

# population_dirs: population folders of the scenario drivers (one per scenario, with gwas/).
# Returns a list of tables, each with the columns H2 and Mc of its scenario:
#   hits         GWAS results with P_Value < p_threshold (SNP, LogP, P_Value, Score, Type)
#   truth        causal SNPs (Shape, Parameter, SNP, Effect)
#   tested       every PPI x shape a GWAS was run for (Score, Type)
#   populations  the folder of every scenario
combine_scenario_runs <- function(population_dirs, p_threshold = 1e-3) {
  parts <- lapply(population_dirs, function(dir) {
    gwas_file <- file.path(dir, "gwas", "gwas_results.csv")
    if (!file.exists(gwas_file)) stop("no gwas/gwas_results.csv in ", dir, " (run the scenario first)")
    cfg <- readRDS(file.path(dir, "data", "config.rds"))
    gw <- read.csv(gwas_file, stringsAsFactors = FALSE)
    tr <- read.csv(file.path(dir, "data", "truth.csv"), stringsAsFactors = FALSE)
    names(tr)[names(tr) == "Trait"] <- "Shape"
    tag <- function(d) cbind(data.frame(H2 = rep(cfg$h2, nrow(d)), Mc = rep(cfg$Mc, nrow(d))), d)
    list(
      hits = tag(gw[gw$P_Value < p_threshold, , drop = FALSE]),
      truth = tag(tr),
      tested = tag(unique(gw[, c("Score", "Type")])),
      populations = data.frame(H2 = cfg$h2, Mc = cfg$Mc, Population = dir, stringsAsFactors = FALSE)
    )
  })
  out <- lapply(c(hits = "hits", truth = "truth", tested = "tested", populations = "populations"),
                function(k) do.call(rbind, lapply(parts, `[[`, k)))
  lapply(out, function(d) { rownames(d) <- NULL; d })
}

# Found causal SNPs per scenario, shape, parameter and PPI, from combine_scenario_runs().
# Each shape is evaluated against its own causal SNPs only (compute_power): a hit on a
# causal SNP of another shape is a false positive. Parameter "False_Positives" holds
# the number of such hits; PPIs without any hit appear with Found = 0.
detection_counts <- function(combined) {
  keys <- unique(combined$tested[, c("H2", "Mc", "Type")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    in_key <- function(d) d$H2 == k$H2 & d$Mc == k$Mc
    hits <- combined$hits[in_key(combined$hits) & combined$hits$Type == k$Type, ]
    truth <- combined$truth[in_key(combined$truth) & combined$truth$Shape == k$Type, ]
    scores <- combined$tested$Score[in_key(combined$tested) & combined$tested$Type == k$Type]
    res <- compute_power(hits, truth, scores)
    cbind(data.frame(H2 = k$H2, Mc = k$Mc, Shape = k$Type), res)
  })
  out <- do.call(rbind, rows)
  out$Found <- as.integer(out$Found)
  rownames(out) <- NULL
  out
}
