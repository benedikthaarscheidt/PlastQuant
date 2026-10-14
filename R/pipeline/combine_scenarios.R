# =============================================================================
# combine_scenarios.R — collect the runs of several scenarios into one table set
# =============================================================================

# Newest run folder (run_<timestamp>) in a scenario's out_dir
latest_run_dir <- function(out_dir) {
  runs <- list.dirs(out_dir, full.names = TRUE, recursive = FALSE)
  runs <- runs[grepl("^run_", basename(runs))]
  if (!length(runs)) stop("no run_* folder in ", out_dir, " (run the scenario first)")
  sort(runs, decreasing = TRUE)[1]
}

# scenario_dirs: out_dirs of run_simulation_scenario(); the newest run of each is used.
# Returns a list of tables, each with the columns H2 and Mc of its scenario:
#   hits    GWAS results with P_Value < p_threshold (SNP, LogP, P_Value, Score, Type)
#   truth   causal SNPs (Shape, Parameter, SNP, Effect)
#   tested  every PPI x shape a GWAS was run for (Score, Type)
#   runs    the run folder of every scenario
combine_scenario_runs <- function(scenario_dirs, p_threshold = 1e-3) {
  parts <- lapply(scenario_dirs, function(dir) {
    run <- latest_run_dir(dir)
    cfg <- readRDS(file.path(run, "config.rds"))
    gw <- read.csv(file.path(run, "gwas_results.csv"), stringsAsFactors = FALSE)
    tr <- read.csv(file.path(run, "truth.csv"), stringsAsFactors = FALSE)
    tag <- function(d) cbind(data.frame(H2 = rep(cfg$h2, nrow(d)), Mc = rep(cfg$Mc, nrow(d))), d)
    list(
      hits = tag(gw[gw$P_Value < p_threshold, , drop = FALSE]),
      truth = tag(tr),
      tested = tag(unique(gw[, c("Score", "Type")])),
      runs = data.frame(H2 = cfg$h2, Mc = cfg$Mc, Run = run, stringsAsFactors = FALSE)
    )
  })
  out <- lapply(c(hits = "hits", truth = "truth", tested = "tested", runs = "runs"),
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
