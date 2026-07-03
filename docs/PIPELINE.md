# PlastQuant pipeline reference

This document describes every script in `R/`, what it does, what it needs, what it
produces, and every parameter you can change. Each script also carries the same
information as a header comment at its top — this file is the cross-script overview.

## Conventions

Every script starts with a header block. Parameters are set with **guarded defaults**:

    if (!exists("SEED")) SEED <- 42

So running a script on its own uses the documented default; running it after a
`scenario_*.R` driver (or another script) that already set the value uses that value.
To change a parameter for a one-off run, edit the value in the script's
`PARAMETERS (edit here)` block, or set it in the R session before sourcing the script.
Always run from the project root (`setwd("~/PlastQuant")` or open the `.Rproj`) so
`here::here()` resolves.

## Dependency map

    01_simulate_genetics ─┐
    02_reaction_norms_*  ─┴─▶ 03_plasticity_scores ─▶ 04_run_gwas ─▶ 05_evaluate_gwas ─▶ 06_analyze_gwas_results
                                     │                                                    └▶ 06_analyze_maize_gwas_results
                                     ├─▶ analysis/summary_stats_vs_scores  (auto-sources 03)
                                     ├─▶ analysis/clustering_rn_summary_stats  (auto-sources 03)
                                     ├─▶ analysis/kendall_rank_corr  (auto-sources 03)
                                     └─▶ analysis/regression_sumstats_dataprep ─▶ analysis/regression_summarystats
    figures/*  read outputs written under output/ by the scripts above.

`03_plasticity_scores.R` auto-sources `01` and `02`, so running it (or anything that
sources it) reruns the upstream simulation. `regression_summarystats.R` requires
`regression_sumstats_dataprep.R` to have been run first (it reads that script's CSV).

## Scripts

<!-- one subsection per script appended by later tasks -->
