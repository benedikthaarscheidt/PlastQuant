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

The simulation model, the PPI computation and the GWAS live as plain functions in
`R/sim`, `R/ppi` and `R/pipeline` (see below); they use no globals. The numbered
scripts and the analyses keep the guarded-default convention.

## Dependency map

    R/sim (simulation model) ─┬─▶ R/pipeline: run_simulation_scenario ◀─ scenarios/scenario_<h2>_<Mc>.R
    R/ppi (PPI registry)     ─┘        (PPIs ─▶ GWAS ─▶ power per shape)
                                              │
                                              ├─▶ 05_evaluate_gwas  (power heatmaps of one run; sourced by the drivers)
                                              └─▶ scenarios/combine_scenarios ─▶ 06_analyze_gwas_results
    R/sim ─▶ 03_plasticity_scores ─▶ 04_run_gwas
                    ├─▶ analysis/summary_stats_vs_scores  (auto-sources 03)
                    ├─▶ analysis/clustering_rn_summary_stats  (auto-sources 03)
                    ├─▶ analysis/kendall_rank_corr  (auto-sources 03)
                    └─▶ analysis/regression_sumstats_dataprep ─▶ analysis/regression_summarystats
    03_plasticity_scores_maize ─▶ 04_run_gwas (scenario_maize) ─▶ 06_analyze_maize_gwas_results
    figures/*  read outputs written under output/ by the scripts above
               (plotting_reaction_norms_all reads a run of scenario_20_5).

The scenario drivers run the new pipeline directly and then 05 for the power heatmaps;
`combine_scenarios.R` collects all scenario runs for 06. The power of a parameter is
always evaluated within its shape: its causal SNPs are searched only in the GWAS of
that shape, and a hit on a SNP that is not causal for the shape is a false positive. `03_plasticity_scores.R` takes its
reaction norms and SNPs from `simulate_reaction_norms()` and keeps the globals that
the analyses and figures read. `regression_summarystats.R` requires
`regression_sumstats_dataprep.R` to have been run first (it reads that script's CSV).

## Scripts

### Simulation model, PPIs and pipeline (R/sim, R/ppi, R/pipeline)

These three folders hold functions only: no global variables, no guarded defaults.
Every function takes its inputs as arguments and returns its result.
`source(here::here("R", "pipeline", "load_pipeline.R"))` loads all of them.

| File | Contents |
|---|---|
| `R/sim/simulation_config.R` | `sim_config()`: every setting, including all variance components; the study space; the parameters of each shape |
| `R/sim/snp_genotypes.R` | SNP dosages, disjoint causal SNP sets, SNP score |
| `R/sim/genetic_values.R` | genetic value of a genotype, repeatability, replicate spread, plant latent values |
| `R/sim/parameter_mapping.R` | latent value → parameter value |
| `R/sim/reaction_norm_shapes.R` | raw shapes, normalisation, reaction-norm curve |
| `R/sim/measurement.R` | measurement error |
| `R/sim/simulate_reaction_norms.R` | `simulate_reaction_norms(cfg)`: runs the model for every shape |
| `R/sim/variance_report.R` | realised variance components on the latent and the trait level |
| `R/ppi/ppi_registry.R` | PPI entries, `register_ppi()`, validation, default registry (28 PPIs) |
| `R/ppi/compute_ppis.R` | replicate means; `compute_ppis()` with `exclude` |
| `R/pipeline/pipeline_input.R` | data contract for simulated and own data |
| `R/pipeline/gwas.R` | kinship matrix; `run_gwas()` per PPI and shape (rrBLUP mixed model) |
| `R/pipeline/evaluate_power.R` | power per shape against that shape's causal SNPs |
| `R/pipeline/run_pipeline.R` | `run_pipeline()`, `run_simulation_scenario()` |
| `R/pipeline/combine_scenarios.R` | `combine_scenario_runs()` (newest run of every scenario), `detection_counts()` per shape |
| `R/pipeline/load_pipeline.R` | sources the three folders |

**Simulation.** `simulate_reaction_norms(sim_config(h2, Mc))` does the following:

1. Simulates the SNP matrix and draws the causal SNPs.
2. Computes, in the order of the model, genetic values, repeatability, plant latent
   values, parameter values, mother reaction norms and replicates, and measurements,
   for each of the four shapes (linear, gaussian, sinusoidal, wave).

All variance settings and their defaults are in `sim_config()`:

- `h2`: share of the latent variance from the causal SNPs.
- `sigma_L`: scale of the level.
- `repeatability_beta`: Beta(18, 2).
- `error_share`: measurement error share, 0.05.

`sim_config()` also holds the remaining settings:

- `n_per_shape` = 200
- `n_rep` = 3
- `n_snps` = 200000
- `n_points` = 6
- `seed`
- `study_space`

The same `seed` gives the same SNPs, causal sets, effects and normal draws for every
`h2`. `variance_report(sim)` shows the realised components:

- On the latent level, the SNP share equals `h2`.
- On the trait level, it shows SNPs, background, interaction, within genotype and
  measurement error.

**Scenarios.** Each `R/scenarios/scenario_<pct>_<Mc>.R` is a single call:

    run_simulation_scenario(sim_config(h2 = 0.2, Mc = 5), out_dir = here::here("output", "scenario_20_5"))

It simulates, computes all PPIs, runs the GWAS with the kinship matrix and evaluates
the power per shape. The first token of the file name is h² in percent; the second
is the number of causal SNPs per parameter.

**Own data.** Build the input and run it:

    input <- pipeline_input(phenotypes, environments, snp_matrix = X, truth = known)
    run_pipeline(input, out_dir)

The input arguments:

- `phenotypes`: a named list with one data.frame per shape or trait. Columns:
  `Genotype`, `Replicate`, `Environment`, `Trait`.
- `environments`: the environment values.
- `snp_matrix`: dosages 0/1/2, with genotype IDs as row names.
- `truth` (optional): known causal SNPs, with columns `Shape`, `Parameter` and `SNP`,
  e.g. known plasticity genes.

Without SNPs, use `run_pipeline(pipeline_input(phenotypes, environments), out_dir, gwas = FALSE)`.
This computes the PPIs only.

**PPIs.** Each PPI is one registry entry that declares its inputs:

- `trait`: trait vector
- `env`: environment values
- `env_index`
- `covariate`
- `matrix`: genotype × environment matrix, at the population level (FW). A
  population-level PPI is always computed once on the replicate means and gives one
  value per genotype, also with `replicate_handling = "replicates"`.

To add a PPI:

    reg <- register_ppi(default_ppi_registry(),
                        ppi_entry("Name", fun, inputs = c(trait_values = "trait", env_values = "env")))

To leave PPIs out, pass `exclude = c("CV_t", ...)`. `replicate_handling` controls
what the GWAS uses:

- `"mean_before"` (default): replicates are averaged before the PPI is computed.
- `"replicates"`: one PPI value per replicate reaction norm, with the line ID repeated.

**Output of a run** (`out_dir/run_<timestamp>/`):

- Always: `phenotypes.rds`, `ppi_scores.rds`, `ppi_scores/<PPI>.csv` (columns `Score`,
  `Type`, `Genotype_ID`, `Replicate`).
- With SNPs: `snp_matrix.rds`.
- With truth: `truth.csv`.
- With GWAS: `gwas_results.csv` (`SNP`, `LogP`, `P_Value`, `Score`, `Type`).
- With truth and GWAS: `gwas_power.csv`, i.e. power per shape, PPI and parameter,
  Bonferroni 0.05 / number of SNPs. Hits on SNPs that are not causal for that shape
  count as false positives.

Scenarios also write `variance_latent.csv`, `variance_trait.csv`, `measurement.csv`,
`mother_parameters.csv`, `mother_curves.rds` and `config.rds`.

**Tests.** `Rscript --vanilla tests/run_tests.R` (from the project root).

### R/03_plasticity_scores.R

```text
=============================================================================
03_plasticity_scores.R — compute every plasticity index (pipeline hub)
=============================================================================
WHAT IT DOES: Simulates the population with the model in R/sim (SNPs -> genetic
  values -> mother norms and replicates -> measurements) and computes every
  plasticity index (via the ppindices package) for each genotype's reaction norm
  across the four forms (linear, gaussian, sinusoidal, wave), assembling the
  `scores_list` consumed by GWAS and analysis.
REQUIRES:     R/sim, R/ppi, R/pipeline (loaded via R/pipeline/load_pipeline.R);
              library(ppindices). The caller MUST set H2, CAUSAL_SNP_NUM and
              OUTPUT_BASE first — the script stops otherwise.
PRODUCES:     `scores_list` in memory; per-score CSVs under OUTPUT_BASE when SAVE=TRUE.
HOW TO RUN:   setwd("~/PlastQuant")
              H2 <- 0.2; CAUSAL_SNP_NUM <- 5
              OUTPUT_BASE <- here::here("output","my_run")
              source(here::here("R", "03_plasticity_scores.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these in a driver script or your R session before sourcing.
  H2                    numeric 0..1, REQUIRED  share of genetic variance from the causal SNPs [COMMON]
  CAUSAL_SNP_NUM        integer, REQUIRED      causal SNPs per parameter                    [COMMON]
  OUTPUT_BASE           path,    REQUIRED      output location                              [COMMON]
  NUM_GENOTYPES         integer, default 800   number of genotypes (divisible by 4)
  NUM_SNPs              integer, default 200000 number of simulated SNPs
  N_MEAS_POINTS         integer, default 6     simulated measurement points per reaction norm
  SEED                  integer, default 1     same seed = same genetics for every H2
  indices, env          optional               subset of the measurement points and their
                                               environment values (default: all points)
  KEEP_REPLICATES       logical, default FALSE keep replicate rows
  GWAS                  logical                run GWAS downstream (read by 04/05/06)
  SAVE                  logical                write per-score CSVs under OUTPUT_BASE
All other model settings (variance components, study space) are the defaults of
sim_config() in R/sim/simulation_config.R.
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
05_evaluate_gwas.R — GWAS accuracy evaluation, per reaction-norm shape
=============================================================================
WHAT IT DOES: Compares the GWAS hits of one run against the simulated causal SNPs.
  Every shape is evaluated against its own causal SNPs only, so the power of a
  parameter is always the power within one shape (a hit on a causal SNP of another
  shape is a false positive). Draws one power heatmap per shape plus an overview.
REQUIRES:     A run folder written by run_simulation_scenario() / run_pipeline()
              (gwas_results.csv, truth.csv, config.rds), e.g. by a scenario_*.R driver.
PRODUCES:     RUN_DIR/gwas_accuracy_summary.csv and RUN_DIR/plots/gwas_power_heatmap.pdf.
HOW TO RUN:   setwd("~/PlastQuant")
              OUTPUT_BASE <- here::here("output", "scenario_20_5")   # newest run is used
              source(here::here("R", "05_evaluate_gwas.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  RUN_DIR      path, default newest run_* in OUTPUT_BASE    run folder to evaluate  [COMMON]
  OUTPUT_BASE  path                                         scenario folder (if no RUN_DIR)
  ALPHA        numeric, default 0.05                        Bonferroni level: P < ALPHA / n_SNPs
=============================================================================
```

### R/06_analyze_gwas_results.R

```text
=============================================================================
06_analyze_gwas_results.R — summarise the GWAS results of all scenarios
=============================================================================
WHAT IT DOES: Summarises, per plasticity index, how many causal SNPs are found in
  every heritability scenario (h2 x causal SNPs per parameter). A parameter is
  always the parameter of one shape ("shape: parameter"); its causal SNPs are
  searched only in the GWAS of that shape, and a hit on a SNP that is not causal
  for that shape is a false positive of that shape. Draws count/fraction heatmaps,
  precision and recall per PPI, true/false positives, and PPI coverage summaries.
REQUIRES:     The tables of R/scenarios/combine_scenarios.R (run the scenarios first).
PRODUCES:     Under COMBINED_DIR: plots/gwas_scenario_summary.pdf,
              precision_recall_per_ppi.csv, recall_per_shape_parameter.csv.
HOW TO RUN:   Rscript --vanilla R/06_analyze_gwas_results.R   (from the repo root)
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  COMBINED_DIR  path, default output/combined_scenarios   tables of combine_scenarios.R [COMMON]
  PPIs          character, default all PPIs in the tables  which plasticity indices to summarise
  FOCUS_PPIS    character, default c("RN", "RTR", "FW")    PPI set for the coverage summary
  OTHER_PPIS    character, default c("APC", "D_slope", "PSI", "RNN")  second PPI set
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
scenario_20_5.R — heritability scenario: h2 = 0.2, 5 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_20_5/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_20_5.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_20_10.R

```text
=============================================================================
scenario_20_10.R — heritability scenario: h2 = 0.2, 10 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_20_10/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_20_10.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_20_20.R

```text
=============================================================================
scenario_20_20.R — heritability scenario: h2 = 0.2, 20 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_20_20/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_20_20.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_40_5.R

```text
=============================================================================
scenario_40_5.R — heritability scenario: h2 = 0.4, 5 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_40_5/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_40_5.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_40_10.R

```text
=============================================================================
scenario_40_10.R — heritability scenario: h2 = 0.4, 10 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_40_10/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_40_10.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_40_20.R

```text
=============================================================================
scenario_40_20.R — heritability scenario: h2 = 0.4, 20 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_40_20/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_40_20.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_60_5.R

```text
=============================================================================
scenario_60_5.R — heritability scenario: h2 = 0.6, 5 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_60_5/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_60_5.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_60_10.R

```text
=============================================================================
scenario_60_10.R — heritability scenario: h2 = 0.6, 10 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_60_10/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_60_10.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_60_20.R

```text
=============================================================================
scenario_60_20.R — heritability scenario: h2 = 0.6, 20 causal SNPs per parameter
=============================================================================
Simulates the population (R/sim), computes all PPIs, runs the GWAS and the
per-shape evaluation, writes everything to output/scenario_60_20/run_<time>/
and draws the power heatmaps (05_evaluate_gwas.R).
All other settings (variance components, study space, ...) are the defaults of
sim_config() in R/sim/simulation_config.R.
HOW TO RUN: Rscript --vanilla R/scenarios/scenario_60_20.R   (from the repo root)
=============================================================================
```

### R/scenarios/scenario_maize.R

```text
=============================================================================
scenario_maize.R — scenario driver: empirical maize dataset
=============================================================================
WHAT IT DOES: Sets the parameters for the empirical maize analysis and runs the GWAS
  on the maize scores (04_run_gwas.R). The maize score computation
  (03_plasticity_scores_maize.R) is available but sourced separately. There is no
  simulated truth for maize, so 05_evaluate_gwas.R does not apply.
REQUIRES:     None to launch — self-contained driver. Uses data/maize_scores_list.rds
              and the maize datasets under data/A-Datasets/.
PRODUCES:     Maize GWAS results under OUTPUT_BASE (output/scenario_maize/).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","scenarios","scenario_maize.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit here — set UNCONDITIONALLY below so they override script defaults)
  OUTPUT_BASE           output location for the maize run
  DO_NOT_PLOT_MAHATTAN  = TRUE  skip Manhattan plots
=============================================================================
```

### R/scenarios/combine_scenarios.R

```text
=============================================================================
combine_scenarios.R — collect the heritability scenarios into one table set
=============================================================================
WHAT IT DOES: Takes the newest run of every scenario_<h2>_<Mc> folder, keeps the
  GWAS hits below P_THRESHOLD and the causal SNPs (both tagged with H2 and Mc), and
  counts the found causal SNPs per scenario, shape, parameter and PPI. Every shape
  is evaluated against its own causal SNPs only (see detection_counts()).
REQUIRES:     Scenario runs written by R/scenarios/scenario_*.R.
PRODUCES:     Under OUT_DIR: gwas_hits.csv, causal_snps_truth.csv, tested.csv,
              detection_counts.csv, runs.csv (read by 06_analyze_gwas_results.R).
HOW TO RUN:   Rscript --vanilla R/scenarios/combine_scenarios.R   (from the repo root)
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  SCENARIO_ROOT  path, default output/                   folder with the scenario_* folders
  SCENARIOS      character, default all scenario_<h2>_<Mc> folders in SCENARIO_ROOT
  P_THRESHOLD    numeric, default 1e-3                   P value below which a SNP is a hit [COMMON]
  OUT_DIR        path, default output/combined_scenarios where the tables are written
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
REQUIRES:     Runs 03 (so it needs H2 / CAUSAL_SNP_NUM / OUTPUT_BASE set
              by a scenario or your session, as 03 does).
PRODUCES:     Regression-input CSV(s) under output/regression_summary_stats/.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","regression_sumstats_dataprep.R"))
-----------------------------------------------------------------------------
PARAMETERS (supplied by 03 / a scenario / your session; see 03_plasticity_scores.R)
  NUM_GENOTYPES, H2, CAUSAL_SNP_NUM, OUTPUT_BASE      [COMMON]
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
REQUIRES:     03 (so H2 / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
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
REQUIRES:     03 (so H2 / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
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
REQUIRES:     03 (so H2 / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
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
WHAT IT DOES: Assembles the reaction-norm overview figure: per shape the mother
  reaction norms of selected genotypes (smooth curve from the simulated mother
  parameters), their replicate measurements (dashed) and replicate range (ribbon),
  plus histograms of the measured trait values.
REQUIRES:     A scenario run (phenotypes.rds, mother_parameters.csv, config.rds),
              by default the newest run of output/scenario_20_5 (run that scenario first).
PRODUCES:     output/plots/reaction_norms_final_figure.pdf,
              reaction_norms_final_figure_panels.pdf, p_combined.pdf.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plotting_reaction_norms_all.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  RUN_DIR         path, default newest run of output/scenario_20_5   input run      [COMMON]
  PLOT_GENOTYPES  list shape -> genotype IDs, default: the genotypes at the 10/50/90 %
                  quantile of the amplitude a of each shape (3 per panel)
  N_CURVE_POINTS  integer, default 200   points of the smooth mother curves
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
REQUIRES:     03 (so H2 / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
PRODUCES:     Resolution-experiment tables/figures under output/.
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS: inherits 03's parameters via the auto-source; resolution knobs are set in
  the body of the script (edit them there).
=============================================================================
```

