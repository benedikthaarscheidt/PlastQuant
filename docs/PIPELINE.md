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

### R/01_simulate_genetics.R

```text
=============================================================================
01_simulate_genetics.R — simulate SNP genotypes and genetic parameters
=============================================================================
WHAT IT DOES: Generates synthetic SNP genotypes for NUM_GENOTYPES individuals with
  optional nested population structure, assigns causal SNPs and effect sizes to the
  reaction-norm parameters, and records the ground-truth causal map later used to
  score GWAS accuracy.
REQUIRES:     Nothing upstream (pure simulation; reads no external data).
PRODUCES:     In memory: the genotype matrix and truth table consumed by 02/03.
              If OUTPUT_BASE is set, writes them under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "01_simulate_genetics.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here)
  NUM_GENOTYPES         integer, default 80      number of simulated genotypes           [COMMON]
  GENETIC_VARIANCES     logical/vector, default FALSE  per-parameter genetic variances
                                                 (FALSE = use built-in defaults)
  POLY_MODE             "uniform"|"exponential", default "uniform"  causal-effect distribution
  POLY_STRENGTH         integer 1-5, default 5   genetic-variance strength (5 = full target SD)
  POLY_COUNTS           integer/vector, default 1  number of causal SNPs per parameter
  STRUCTURED_POPULATION logical, default FALSE   add nested population structure
  SEED                  integer, default 42      RNG seed for reproducibility            [COMMON]
  OUTPUT_BASE           path, default here::here("output","default")  where outputs go
-----------------------------------------------------------------------------
These parameters are supplied by the sourcing script (03_plasticity_scores.R or a
scenario_*.R driver) or by your R session before sourcing. They are intentionally
NOT assigned here: the simulation branches on whether SEED / GENETIC_VARIANCES exist
(see below), so defining them here would change the RNG path. GENETIC_VARIANCES falls
back to FALSE further down; SEED, when defined, offsets the seeds (leave it unset to
reproduce the default simulation).
=============================================================================
```

### R/02_reaction_norms_linear.R

```text
=============================================================================
02_reaction_norms_linear.R — generate linear reaction norms
=============================================================================
WHAT IT DOES: Builds linear (intercept + slope) reaction norms across environments
  for each genotype, optionally driven by the simulated genetics from 01.
REQUIRES:     Genetic parameters from 01_simulate_genetics.R when USE_GENETICS = TRUE
              (typically sourced upstream by 03_plasticity_scores.R).
PRODUCES:     In memory: the linear reaction-norm data consumed by 03.
HOW TO RUN:   Sourced by 03_plasticity_scores.R; or, with 01 already run,
              setwd("~/PlastQuant"); source(here::here("R", "02_reaction_norms_linear.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here)
  USE_GENETICS  logical, default FALSE   use simulated genetics (FALSE = random params)  [COMMON]
  SEED          integer, default 42      RNG seed                                        [COMMON]
  OUTPUT_BASE   path, default here::here("output","default")  where outputs go
-----------------------------------------------------------------------------
These parameters are supplied by the sourcing script (03_plasticity_scores.R) or your
R session. They are intentionally NOT assigned here: the script branches on whether
SEED / USE_GENETICS exist (USE_GENETICS falls back to FALSE further down; SEED, when
defined, offsets the RNG seed), so defining them here would change behaviour.
=============================================================================
```

### R/02_reaction_norms_nonlinear.R

```text
=============================================================================
02_reaction_norms_nonlinear.R — generate nonlinear reaction norms
=============================================================================
WHAT IT DOES: Builds nonlinear reaction norms (Gaussian, sinusoidal, and wave
  shapes) across environments for each genotype, optionally driven by the
  simulated genetics from 01.
REQUIRES:     Genetic parameters from 01_simulate_genetics.R when USE_GENETICS = TRUE
              (typically sourced upstream by 03_plasticity_scores.R).
PRODUCES:     In memory: the nonlinear reaction-norm data consumed by 03.
HOW TO RUN:   Sourced by 03_plasticity_scores.R; or, with 01 already run,
              setwd("~/PlastQuant"); source(here::here("R", "02_reaction_norms_nonlinear.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here)
  USE_GENETICS  logical, default FALSE   use simulated genetics (FALSE = random params)  [COMMON]
  SEED          integer, default 42      RNG seed                                        [COMMON]
  OUTPUT_BASE   path, default here::here("output","default")  where outputs go
-----------------------------------------------------------------------------
These parameters are supplied by the sourcing script (03_plasticity_scores.R) or your
R session. They are intentionally NOT assigned here: the script branches on whether
SEED / USE_GENETICS exist (USE_GENETICS falls back to FALSE further down; SEED, when
defined, offsets the RNG seed), so defining them here would change behaviour.
=============================================================================
```

### R/03_plasticity_scores.R

```text
=============================================================================
03_plasticity_scores.R — compute every plasticity index (pipeline hub)
=============================================================================
WHAT IT DOES: Computes every plasticity index (via the ppindices package) for each
  simulated genotype's reaction norm across the four forms (linear, gaussian,
  sinusoidal, wave), assembling the `scores_list` consumed by GWAS and analysis.
  This is the hub: it auto-sources 01 and 02, so running it reruns the upstream
  simulation.
REQUIRES:     Auto-sources 01_simulate_genetics.R, 02_reaction_norms_linear.R,
              02_reaction_norms_nonlinear.R; library(ppindices).
              When USE_GENETICS = TRUE (the default) the caller MUST set
              HERITABILITY_5TH, CAUSAL_SNP_NUM and OUTPUT_BASE first (a scenario_*.R
              driver, or assignments in your R session) — the script stops otherwise.
PRODUCES:     `scores_list` in memory; per-score CSVs under OUTPUT_BASE when SAVE=TRUE.
HOW TO RUN:   setwd("~/PlastQuant")
              HERITABILITY_5TH <- 1; CAUSAL_SNP_NUM <- 5
              OUTPUT_BASE <- here::here("output","my_run")
              source(here::here("R", "03_plasticity_scores.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these in a scenario_*.R driver or your R session before sourcing.
  The defaults below are applied in-script (guarded) at the noted line numbers.
  USE_GENETICS          logical, default TRUE  use simulated genetics (line 36)        [COMMON]
  HERITABILITY_5TH      numeric, REQUIRED when USE_GENETICS=TRUE  H^2 of the 5th trait [COMMON]
  CAUSAL_SNP_NUM        integer, REQUIRED when USE_GENETICS=TRUE  causal SNPs / param  [COMMON]
  OUTPUT_BASE           path,    REQUIRED when USE_GENETICS=TRUE  output location      [COMMON]
  NUM_GENOTYPES         integer, default 800   number of genotypes (line 100)
  N_MEAS_POINTS         integer, default 6     measurement points per reaction norm (line 161)
  KEEP_REPLICATES       logical, default (see lines 97/179)  keep replicate rows
  STRUCTURED_POPULATION logical, default FALSE nested population structure (line 98)
  GENETIC_VARIANCES     logical/vector         per-parameter genetic variances (from 01)
  POLY_MODE/STRENGTH/COUNTS  causal-effect controls passed through to 01
  GWAS                  logical                run GWAS downstream (read by 04/05/06)
  SAVE                  logical                write per-score CSVs under OUTPUT_BASE
=============================================================================
```

### R/03_plasticity_scores_maize.R

```text
=============================================================================
03_plasticity_scores_maize.R — compute plasticity indices for the maize data
=============================================================================
WHAT IT DOES: Loads the empirical maize datasets and computes every plasticity index
  (via the ppindices package) across the three phenotypes (biomass, leaf area, water-
  use efficiency), producing the maize `scores_list`. Maize counterpart of
  03_plasticity_scores.R (loads real data instead of simulating it).
REQUIRES:     data/A-Datasets/ maize inputs; library(ppindices).
PRODUCES:     Maize `scores_list` in memory; per-score CSVs under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "03_plasticity_scores_maize.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top / set in your session)
  OUTPUT_BASE           path, default output/scenario_maize  where outputs go; edit `OUTPUT_BASE <-`
  NUM_GENOTYPES         integer, default 800  guarded default (from data if provided)
  NUM_SNPs              integer, default 200000  guarded default
  USE_GENETICS          logical, default TRUE  guarded further down
  KEEP_REPLICATES / STRUCTURED_POPULATION / GENETIC_VARIANCES  guarded defaults below
=============================================================================
```

### R/04_run_gwas.R

```text
=============================================================================
04_run_gwas.R — GWAS on plasticity scores
=============================================================================
WHAT IT DOES: Runs GWAS (rrBLUP / sommer) for each plasticity score against the
  simulated genotypes, optionally with a kinship correction, and (unless suppressed)
  writes Manhattan plots. Identifies which SNPs each index recovers.
REQUIRES:     scores_list and genetics from 03_plasticity_scores.R (run 03 first, or
              source it); data/A-Datasets/ genotype inputs; packages rrBLUP / sommer.
PRODUCES:     GWAS result CSVs and Manhattan PDFs under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "04_run_gwas.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit at the noted line)
  USE_KINSHIP           logical, default TRUE   include kinship (K) in the mixed model;
                                                edit the `USE_KINSHIP <-` assignment below [COMMON]
  DO_NOT_PLOT_MAHATTAN  if defined, Manhattan PDFs are skipped (faster)                 [COMMON]
  GWAS                  logical (from caller/03) whether GWAS is run                     [COMMON]
  OUTPUT_BASE           path (from caller/03)    where results are written
=============================================================================
```

### R/05_evaluate_gwas.R

```text
=============================================================================
05_evaluate_gwas.R — GWAS accuracy evaluation
=============================================================================
WHAT IT DOES: Compares GWAS hits against the simulated ground-truth causal SNPs to
  quantify how well each plasticity index recovers the true signal. Produces 5 power
  heatmaps: one combined (all forms) plus one per reaction-norm form.
REQUIRES:     GWAS results from 04_run_gwas.R — this script auto-sources 04 (which in
              turn needs 03). Reads the GWAS results CSV under OUTPUT_BASE.
PRODUCES:     5 power-heatmap PDFs under OUTPUT_BASE/plots.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "05_evaluate_gwas.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit at the noted line)
  alpha   numeric, default 0.05   significance threshold; edit the `alpha <-` assignment  [COMMON]
  GWAS    logical (from caller)    whether GWAS was run upstream
=============================================================================
```

### R/06_analyze_gwas_results.R

```text
=============================================================================
06_analyze_gwas_results.R — summarise simulated-GWAS results
=============================================================================
WHAT IT DOES: Loads the causal-SNP ground truth and the simulated GWAS result tables
  and summarises, per plasticity index, how the recovered signal compares to truth,
  rendering summary heatmaps/plots.
REQUIRES:     GWAS result CSVs produced by 04/05 (run a scenario first); the causal
              truth text files under the scenario's OUTPUT_BASE.
PRODUCES:     Summary heatmaps/plots under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "06_analyze_gwas_results.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit at the noted line)
  PPIs   character vector   which plasticity indices to summarise; edit `PPIs <-`     [COMMON]
  GWAS   logical            whether GWAS results are present
=============================================================================
```

### R/06_analyze_maize_gwas_results.R

```text
=============================================================================
06_analyze_maize_gwas_results.R — summarise empirical maize-GWAS results
=============================================================================
WHAT IT DOES: Summarises the GWAS results for the empirical maize dataset across the
  selected plasticity indices (PPIs) and reaction-norm traits (RNs), rendering the
  maize result heatmaps/plots.
REQUIRES:     Maize GWAS result data (produced by the maize scenario / 03_maize + 04);
              data/maize_scores_list.rds.
PRODUCES:     Maize summary heatmaps/plots under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R", "06_analyze_maize_gwas_results.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit at the noted line)
  PPIs   character vector   plasticity indices to summarise; edit `PPIs <-`            [COMMON]
  RNs    character vector   reaction-norm traits (leafArea/WUE/biomass); edit `RNs <-` [COMMON]
  GWAS   logical            whether GWAS results are present
=============================================================================
```

### R/scenarios/scenario_20_5.R

```text
=============================================================================
scenario_20_5.R — scenario driver: HERITABILITY_5TH=1, CAUSAL_SNP_NUM=5
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_20_5/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_20_5.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 1    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 5    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_20_10.R

```text
=============================================================================
scenario_20_10.R — scenario driver: HERITABILITY_5TH=1, CAUSAL_SNP_NUM=10
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_20_10/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_20_10.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 1    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 10    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_20_20.R

```text
=============================================================================
scenario_20_20.R — scenario driver: HERITABILITY_5TH=1, CAUSAL_SNP_NUM=20
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_20_20/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_20_20.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 1    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 20    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_40_5.R

```text
=============================================================================
scenario_40_5.R — scenario driver: HERITABILITY_5TH=2, CAUSAL_SNP_NUM=5
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_40_5/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_40_5.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 2    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 5    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_40_10.R

```text
=============================================================================
scenario_40_10.R — scenario driver: HERITABILITY_5TH=2, CAUSAL_SNP_NUM=10
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_40_10/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_40_10.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 2    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 10    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_40_20.R

```text
=============================================================================
scenario_40_20.R — scenario driver: HERITABILITY_5TH=2, CAUSAL_SNP_NUM=20
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_40_20/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_40_20.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 2    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 20    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_60_5.R

```text
=============================================================================
scenario_60_5.R — scenario driver: HERITABILITY_5TH=3, CAUSAL_SNP_NUM=5
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_60_5/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_60_5.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 3    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 5    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_60_10.R

```text
=============================================================================
scenario_60_10.R — scenario driver: HERITABILITY_5TH=3, CAUSAL_SNP_NUM=10
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_60_10/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_60_10.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 3    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 10    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_60_20.R

```text
=============================================================================
scenario_60_20.R — scenario driver: HERITABILITY_5TH=3, CAUSAL_SNP_NUM=20
=============================================================================
WHAT IT DOES: Sets the parameters for this scenario, then runs the plasticity-score
  computation (03_plasticity_scores.R) and the GWAS accuracy evaluation
  (05_evaluate_gwas.R) end-to-end.
REQUIRES:     None — self-contained driver. Run from the project root.
PRODUCES:     Scores + GWAS evaluation under OUTPUT_BASE (output/scenario_60_20/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_60_20.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for this scenario
  HERITABILITY_5TH      = 3    heritability setting for the 5th trait          [COMMON]
  CAUSAL_SNP_NUM        = 20    causal SNPs per parameter                       [COMMON]
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/scenario_maize.R

```text
=============================================================================
scenario_maize.R — scenario driver: empirical maize dataset
=============================================================================
WHAT IT DOES: Sets the parameters for the empirical maize analysis and runs the GWAS
  accuracy evaluation (05_evaluate_gwas.R). The maize score computation
  (03_plasticity_scores_maize.R) is available but sourced separately.
REQUIRES:     None to launch — self-contained driver. Uses data/maize_scores_list.rds
              and the maize datasets under data/A-Datasets/.
PRODUCES:     Maize GWAS evaluation under OUTPUT_BASE (output/scenario_maize/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_maize.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for the maize run
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/analysis/regression_sumstats_dataprep.R

```text
=============================================================================
regression_sumstats_dataprep.R — build the summary-stats vs scores table
=============================================================================
WHAT IT DOES: Assembles the table of reaction-norm summary statistics against
  plasticity scores that regression_summarystats.R and several figure scripts consume.
  Auto-sources 03_plasticity_scores.R and analysis/summary_stats_vs_scores.R.
REQUIRES:     Runs 03 (so it needs HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE set
              by a scenario or your session, as 03 does).
PRODUCES:     Regression-input CSV(s) under output/regression_summary_stats/.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","regression_sumstats_dataprep.R"))
-----------------------------------------------------------------------------
PARAMETERS (supplied by 03 / a scenario / your session; see 03_plasticity_scores.R)
  NUM_GENOTYPES, HERITABILITY_5TH, CAUSAL_SNP_NUM, USE_GENETICS, OUTPUT_BASE      [COMMON]
=============================================================================
```

### R/analysis/regression_summarystats.R

```text
=============================================================================
regression_summarystats.R — regress plasticity scores on summary statistics
=============================================================================
WHAT IT DOES: For each plasticity score, fits full/simple/ridge regression models on
  the nine summary-stat predictors, builds coefficient tables with significance stars
  and R^2, clusters scores, and writes per-score and combined diagnostic PDFs.
REQUIRES:     regression_sumstats_dataprep.R MUST be run first — this script reads its
              combined CSV from output/regression_summary_stats/.
PRODUCES:     Coefficient tables and diagnostic PDFs under output/plots/.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","regression_summarystats.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  input_folder     path, default output/regression_summary_stats  where the CSV is read [COMMON]
  output_folder    path, default output/plots                     where PDFs are written [COMMON]
  predictor_names  character vector  the summary-stat predictors used
  k                integer, default 4     number of score clusters
  n_perm           integer, default 1000  permutations for the cluster test
=============================================================================
```

### R/analysis/kendall_rank_corr.R

```text
=============================================================================
kendall_rank_corr.R — Kendall rank correlations between indices
=============================================================================
WHAT IT DOES: Computes Kendall rank correlations between the plasticity indices over
  the chosen sampling intervals and reaction-norm ranges. Auto-sources
  03_plasticity_scores.R when the scores are not already loaded.
REQUIRES:     03 (so HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
PRODUCES:     Rank-correlation tables (and inputs for downstream figures).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","kendall_rank_corr.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  sampling_intervals    integer vector, default c(10)   measurement intervals to test   [COMMON]
  global_initial_length integer, default 50             full reaction-norm length
  range_list            list, default list(full=c(1,50)) sub-ranges to evaluate
  score_names           character vector                which indices to correlate
=============================================================================
```

### R/analysis/kendall_tau_corr_clustering.R

```text
=============================================================================
kendall_tau_corr_clustering.R — cluster indices by rank agreement
=============================================================================
WHAT IT DOES: Reads the rank-agreement ("rank flops") tables and clusters the
  plasticity indices by their Kendall-tau agreement across top-k rankings.
REQUIRES:     Rank-flop CSVs under output/rank_flops/ (produced upstream).
PRODUCES:     Clustering summaries/plots of index agreement.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","kendall_tau_corr_clustering.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  in_dir          path, default output/rank_flops   input rank-agreement files          [COMMON]
  tau_threshold   numeric, default 0.8              Kendall-tau agreement threshold      [COMMON]
  topk_values     integer vector, default c(3,5,10,15,20)  top-k cutoffs to evaluate     [COMMON]
  AGREE_THR_FULL  numeric, default 1.00             full-agreement threshold
=============================================================================
```

### R/analysis/clustering_rn_summary_stats.R

```text
=============================================================================
clustering_rn_summary_stats.R — cluster reaction-norm summary statistics
=============================================================================
WHAT IT DOES: Clusters genotypes / indices on their reaction-norm summary statistics
  (k-means, DBSCAN, hierarchical, mclust) and renders the associated heatmaps/plots.
  Auto-sources 03_plasticity_scores.R when scores are not already loaded
  (guarded by the `data_loaded` flag at the top).
REQUIRES:     03 (so HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
PRODUCES:     Clustering figures under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","clustering_rn_summary_stats.R"))
-----------------------------------------------------------------------------
PARAMETERS: inherits 03's parameters via the auto-source; clustering knobs are set
  in the body of the script (edit them there). Set `data_loaded <- TRUE` beforehand to
  reuse already-computed scores instead of re-running 03.
=============================================================================
```

### R/analysis/summary_stats_vs_scores.R

```text
=============================================================================
summary_stats_vs_scores.R — relate summary statistics to plasticity scores
=============================================================================
WHAT IT DOES: Computes and visualises the relationship between reaction-norm summary
  statistics and the plasticity scores (correlation matrices, distance measures,
  clustering). Auto-sources 03_plasticity_scores.R when scores are not already loaded
  (guarded by the `data_loaded` flag at the top).
REQUIRES:     03 (so HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
PRODUCES:     Summary-stat vs score figures/tables under OUTPUT_BASE.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","summary_stats_vs_scores.R"))
-----------------------------------------------------------------------------
PARAMETERS: inherits 03's parameters via the auto-source; analysis knobs are set in the
  body. Set `data_loaded <- TRUE` beforehand to reuse already-computed scores.
=============================================================================
```

### R/figures/combined_resolution_figure.R

```text
=============================================================================
combined_resolution_figure.R — combined resolution figure
=============================================================================
WHAT IT DOES: Combines the resolution heatmap (left) and the confounder/summary-stat
  panel (right) into one publication figure.
REQUIRES:     Regression-summary CSVs under output/regression_summary_stats/
              (produced by regression_sumstats_dataprep.R / regression_summarystats.R).
PRODUCES:     A combined resolution figure PDF.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","combined_resolution_figure.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  directory   path, default output/regression_summary_stats  input CSV folder           [COMMON]
  files       character vector  the specific input CSVs to combine
  n_samples   integer vector, default c(2,3,4,5,6,11,26,50)   sample counts per resolution
  stat_cols   character vector  summary-stat columns to use
=============================================================================
```

### R/figures/correlation_resolution.R

```text
=============================================================================
correlation_resolution.R — correlation vs sampling resolution
=============================================================================
WHAT IT DOES: Plots how index correlations change with the sampling resolution of the
  reaction norm.
REQUIRES:     Regression-summary CSVs under output/regression_summary_stats/.
PRODUCES:     Correlation-vs-resolution figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","correlation_resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  directory    path, default output/regression_summary_stats  input CSV folder          [COMMON]
  resolutions  integer vector, default c(50,25,20,15,10,5,2,1) sampling resolutions
  n_samples    integer vector, default c(2,3,4,5,6,11,26,50)   sample counts
  sel_idx      integer vector, default c(1,3,7,8)   which resolutions to plot            [COMMON]
  files        character vector  the input CSVs
  stat_cols    character vector  summary-stat columns to use
=============================================================================
```

### R/figures/correlation_resolution_experiment.R

```text
=============================================================================
correlation_resolution_experiment.R — resolution robustness with confounders
=============================================================================
WHAT IT DOES: A resolution-robustness analysis like correlation_resolution.R, but with
  the summary statistics included in the regression as confounders.
REQUIRES:     Regression-summary CSVs under output/regression_summary_stats/.
PRODUCES:     Resolution-experiment figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","correlation_resolution_experiment.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  directory   path, default output/regression_summary_stats  input CSV folder           [COMMON]
  files       character vector  the input CSVs
  n_samples   integer vector, default c(2,3,4,5,6,11,26,50)   sample counts
  stat_cols   character vector  summary-stat columns to use
=============================================================================
```

### R/figures/correlation_summary_stats_PS.R

```text
=============================================================================
correlation_summary_stats_PS.R — summary-stat vs plasticity-score correlations
=============================================================================
WHAT IT DOES: Renders correlation panels between the reaction-norm summary statistics
  and the plasticity scores.
REQUIRES:     A regression-data CSV under output/regression_summary_stats/.
PRODUCES:     Correlation-panel figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","correlation_summary_stats_PS.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  the input CSV path (read via read_csv at the top)                                      [COMMON]
  stat_cols   character vector  summary-stat columns to use
=============================================================================
```

### R/figures/plasticity_score_realisation_plot.R

```text
=============================================================================
plasticity_score_realisation_plot.R — plasticity-score realisations (simulated)
=============================================================================
WHAT IT DOES: Plots the realised plasticity scores across genotypes/forms for the
  simulated data.
REQUIRES:     Score outputs produced by 03_plasticity_scores.R (run a scenario first),
              read from under output/.
PRODUCES:     Plasticity-score realisation figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plasticity_score_realisation_plot.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named input path / selection assignments in the body)
  the input score path and any index/form selection are set in the script body.
=============================================================================
```

### R/figures/plasticity_score_realisation_plot_maize.R

```text
=============================================================================
plasticity_score_realisation_plot_maize.R — plasticity-score realisations (maize)
=============================================================================
WHAT IT DOES: Plots the realised plasticity scores for the empirical maize dataset.
  Sources figures/rank_corr_biplot.R.
REQUIRES:     Maize score outputs (data/maize_scores_list.rds and the maize pipeline
              outputs); figures/rank_corr_biplot.R.
PRODUCES:     Maize plasticity-score realisation figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plasticity_score_realisation_plot_maize.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  form_ranges  list  the reaction-norm ranges/forms plotted                              [COMMON]
  the input score path is set in the script body.
=============================================================================
```

### R/figures/plotting_reaction_norms_all.R

```text
=============================================================================
plotting_reaction_norms_all.R — reaction-norm overview figure
=============================================================================
WHAT IT DOES: Assembles the reaction-norm overview figure: genotype line plots plus
  histograms of the trait distribution.
REQUIRES:     all_combined_data.csv written by a scenario run under
              output/<scenario>/synthetic_data/... (a scenario must be run first).
PRODUCES:     Reaction-norm overview figure PDF.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plotting_reaction_norms_all.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  the input all_combined_data path (read via read_csv at the top)                        [COMMON]
  ylim   numeric length-2  y-axis limits for the line plots
=============================================================================
```

### R/figures/plotting_resolution.R

```text
=============================================================================
plotting_resolution.R — resolution correlation plots
=============================================================================
WHAT IT DOES: Reads the per-resolution correlation CSVs and renders the resolution
  correlation plots across the defined reaction-norm ranges.
REQUIRES:     correlations_*.csv files under the input directory (produced upstream).
PRODUCES:     Resolution correlation figure(s).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plotting_resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  dir_path     path   folder holding the correlations_*.csv inputs                       [COMMON]
  range_names  character vector, default c("full","partial1","partial2")  ranges plotted
=============================================================================
```

### R/figures/rank_corr_biplot.R

```text
=============================================================================
rank_corr_biplot.R — rank-correlation biplot
=============================================================================
WHAT IT DOES: Builds the rank-correlation biplot of the plasticity indices (also used
  as a helper by plasticity_score_realisation_plot_maize.R).
REQUIRES:     Rank-correlation inputs produced upstream (read from under output/).
PRODUCES:     Rank-correlation biplot figure.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","rank_corr_biplot.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named input path / selection assignments in the body)
  the input path and index selection are set in the script body.
=============================================================================
```

### R/figures/resolution.R

```text
=============================================================================
resolution.R — resolution experiment (dataprep + regression combined)
=============================================================================
WHAT IT DOES: Runs, in one pass, what regression_sumstats_dataprep.R and
  regression_summarystats.R do together, for the resolution experiment. Auto-sources
  03_plasticity_scores.R and analysis/summary_stats_vs_scores.R.
REQUIRES:     03 (so HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
PRODUCES:     Resolution-experiment tables/figures under output/.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS: inherits 03's parameters via the auto-source; resolution knobs are set in
  the body of the script (edit them there).
=============================================================================
```

