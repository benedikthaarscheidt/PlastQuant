# PlastQuant pipeline reference

This document describes every script in `R/`, what it does, what it needs, what it
produces, and every parameter you can change. Each script also carries the same
information as a header comment at its top — this file is the cross-script overview.
For the procedure as a whole, start with the README.

## Conventions

Always run from the project root (`setwd("~/PlastQuant")` or open the `.Rproj`) so
`here::here()` resolves.

The simulation model, the PPI computation, the GWAS and the folder handling are plain
functions in `R/sim`, `R/ppi` and `R/pipeline`; they use no globals.
`source(here::here("R", "pipeline", "load_pipeline.R"))` loads all of them.

The scripts work on folders, set before sourcing:

- `POPULATION_DIR`: a population folder `output/populations/<name>/` (simulation, the
  analyses across resolutions);
- `RESOLUTION_DIR`: one resolution of a population, `output/populations/<name>/n_<k>/`
  (the figures and analyses of one measurement, the GWAS evaluation).

Other parameters of the numbered scripts and analyses are set with guarded defaults
(`if (!exists("X")) X <- ...`) or edited in the script's PARAMETERS block.

## Population and resolution folders

    output/populations/<name>/
      data/          the population: norms.rds (true mother and replicate norms),
                     mother_parameters.csv, plant_parameters.csv, reaction_norms.pdf,
                     snp_matrix.rds|.csv, kinship.csv, truth.csv, config.rds
      variance/      variance_latent.csv
      n_<k>/         the population measured at k equidistant environments
        data/        traits/<shape>.csv (Genotype, one column per environment), environments.csv
        variance/    variance_trait.csv, measurement.csv
        ppis/        <PPI>.csv, ppi_scores.rds
        gwas/        gwas_results.csv, gwas_power.csv, gwas_accuracy_summary.csv
        analyses/    summary_statistics.csv, analyses of this resolution
        plots/       figures of this resolution
      analyses/      analyses across resolutions (resolution/, correlation_summary_stats/)
      plots/         figures across resolutions

Own data: SNPs, kinship and truth in `<name>/data/`, the measurements in
`<name>/n_<number of environments>/data/` (same format as above).

## Dependency map

    simulate_population(cfg, POPULATION_DIR)          -> data/, variance/
      └─ sample_resolution(POPULATION_DIR, k)         -> n_<k>/data/, n_<k>/variance/
           └─ compute_population_ppis(n_<k>)          -> n_<k>/ppis/
                ├─ run_population_gwas(n_<k>)         -> n_<k>/gwas/ ─▶ 05_evaluate_gwas
                ├─ figures/plasticity_score_realisation_plot, rank_corr_biplot,
                │  correlation_summary_stats_PS, plotting_reaction_norms_all
                ├─ analysis/regression_summarystats, kendall_rank_corr ─▶ kendall_tau_corr_clustering
                └─ 03_plasticity_scores ─▶ analysis/summary_stats_vs_scores, clustering_rn_summary_stats

    analysis/regression_sumstats_dataprep (POPULATION_DIR): n_2 ... n_50 (the three steps above per k)
      ├─ analysis/resolution_analysis ─▶ figures/combined_resolution_figure
      └─ figures/correlation_resolution, correlation_resolution_experiment
    figures/resolution (POPULATION_DIR, samples the gradient halves itself) ─▶ figures/plotting_resolution

    scenarios/scenario_<h2>_<Mc>: simulate -> n_6 -> PPIs -> GWAS -> 05
      └─ scenarios/combine_scenarios ─▶ 06_analyze_gwas_results
    03_plasticity_scores_maize ─▶ 04_run_gwas (scenario_maize) ─▶ 06_analyze_maize_gwas_results

The power of a parameter is always evaluated within its shape: its causal SNPs are
searched only in the GWAS of that shape, and a hit on a SNP that is not causal for the
shape is a false positive.

## Functions (R/sim, R/ppi, R/pipeline)

| File | Contents |
|---|---|
| `R/sim/simulation_config.R` | `sim_config()`: every setting, including all variance components (each commented); the study space; the parameters of each shape |
| `R/sim/snp_genotypes.R` | SNP dosages, disjoint causal SNP sets, SNP score |
| `R/sim/genetic_values.R` | genetic value of a genotype, repeatability, replicate spread, plant latent values |
| `R/sim/parameter_mapping.R` | latent value → parameter value |
| `R/sim/reaction_norm_shapes.R` | raw shapes, normalisation, reaction-norm curve |
| `R/sim/measurement.R` | measurement error |
| `R/sim/simulate_reaction_norms.R` | `simulate_norms()` (true norms), `measure_norms()` (sampling at n points or at points t, with measurement error), `simulate_reaction_norms()` (both) |
| `R/sim/variance_report.R` | `latent_variance()` (from the norms), `trait_variance()` (per measurement), `variance_report()` |
| `R/ppi/ppi_registry.R` | PPI entries, `register_ppi()`, validation, default registry (28 PPIs) |
| `R/ppi/compute_ppis.R` | replicate means; `compute_ppis()` with `exclude` and `replicate_handling`; `ppi_wide_table()` |
| `R/pipeline/pipeline_input.R` | in-memory data contract (phenotypes, environments, SNPs, truth) and its checks |
| `R/pipeline/population_folder.R` | folder format: `write_population()`, `write_resolution()`, `read_population()`, `sample_population()`, `resolution_dirs()`, `resolution_path()` |
| `R/pipeline/population_steps.R` | the steps: `simulate_population()`, `sample_resolution()`, `compute_population_ppis()`, `run_population_gwas()` |
| `R/pipeline/plot_norms.R` | `plot_true_norms()`: the true norms, drawn into `data/reaction_norms.pdf` |
| `R/pipeline/gwas.R` | kinship matrix; `run_gwas()` per PPI and shape (rrBLUP mixed model) |
| `R/pipeline/evaluate_power.R` | power per shape against that shape's causal SNPs |
| `R/pipeline/resolution_analysis.R` | summary statistics; `compute_resolutions()`, `resolution_score_tables()`, `population_score_table()`; correlations per resolution, their regression on n, the mixed model |
| `R/pipeline/similarity.R` | similarity measures and their permutation test |
| `R/pipeline/combine_scenarios.R` | `combine_scenario_runs()` (the resolution n_6 of every scenario), `detection_counts()` per shape |
| `R/pipeline/load_pipeline.R` | sources the three folders |

**PPIs.** Each PPI is one registry entry that declares its inputs (`trait`, `env`,
`env_index`, `covariate`, or `matrix` = genotype × environment matrix for the
population-level FW). A population-level PPI is always computed once on the replicate
means and gives one value per genotype. To add a PPI:

    reg <- register_ppi(default_ppi_registry(),
                        ppi_entry("Name", fun, inputs = c(trait_values = "trait", env_values = "env")))
    compute_population_ppis(RESOLUTION_DIR, registry = reg)

**Tests.** `Rscript --vanilla tests/run_tests.R` (from the project root).

## Scripts

### R/03_plasticity_scores.R

```text
=============================================================================
03_plasticity_scores.R — compute every plasticity index (pipeline hub)
=============================================================================
WHAT IT DOES: Loads a resolution of a simulated population (RESOLUTION_DIR, made with
  simulate_population() and sample_resolution()) and computes every plasticity index
  (via the ppindices package) for each genotype's reaction norm across the four forms
  (linear, gaussian, sinusoidal, wave), assembling the `scores_list` consumed by GWAS
  and analysis. Everything it writes goes to RESOLUTION_DIR/analyses/scores/.
REQUIRES:     R/sim, R/ppi, R/pipeline (loaded via R/pipeline/load_pipeline.R);
              library(ppindices). The caller MUST set RESOLUTION_DIR.
PRODUCES:     `scores_list` in memory; outputs under RESOLUTION_DIR/analyses/scores/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "my_population", "n_6")
              source(here::here("R", "03_plasticity_scores.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these in a driver script or your R session before sourcing.
  RESOLUTION_DIR        path,    REQUIRED      resolution folder of a simulated population  [COMMON]
  indices, env          optional               subset of the measurement points and their
                                               environment values (default: all points)
  KEEP_REPLICATES       logical, default FALSE keep replicate rows
  GWAS                  logical                run GWAS downstream (read by 04/05/06)
  SAVE                  logical                write per-score CSVs under analyses/scores/
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
WHAT IT DOES: Compares the GWAS hits of one resolution of a population against the
  population's causal SNPs.
  Every shape is evaluated against its own causal SNPs only, so the power of a
  parameter is always the power within one shape (a hit on a causal SNP of another
  shape is a false positive). Draws one power heatmap per shape plus an overview.
REQUIRES:     A resolution folder with GWAS results (gwas/gwas_results.csv) of a
              population with a truth table (<population>/data/truth.csv), e.g. from a scenario_*.R driver
              or run_population_gwas() on own data with known causal SNPs.
PRODUCES:     RESOLUTION_DIR/gwas/gwas_accuracy_summary.csv and
              RESOLUTION_DIR/plots/gwas_power_heatmap.pdf.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "scenario_20_5", "n_6")
              source(here::here("R", "05_evaluate_gwas.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  RESOLUTION_DIR  path, REQUIRED              resolution folder to evaluate         [COMMON]
  ALPHA           numeric, default 0.05       Bonferroni level: P < ALPHA / n_SNPs
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
Simulates the population (R/sim) into output/populations/scenario_20_5/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_20_10/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_20_20/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_40_5/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_40_10/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_40_20/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_60_5/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_60_10/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
Simulates the population (R/sim) into output/populations/scenario_60_20/ (data/,
variance/), measures it at 6 points (n_6/), computes all PPIs (n_6/ppis/), runs the
GWAS and the per-shape evaluation (n_6/gwas/) and draws the power heatmaps
(n_6/plots/, 05_evaluate_gwas.R).
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
WHAT IT DOES: Takes the population folder of every scenario_<h2>_<Mc>, keeps the
  GWAS hits below P_THRESHOLD and the causal SNPs (both tagged with H2 and Mc), and
  counts the found causal SNPs per scenario, shape, parameter and PPI. Every shape
  is evaluated against its own causal SNPs only (see detection_counts()).
REQUIRES:     Scenario runs written by R/scenarios/scenario_*.R.
PRODUCES:     Under OUT_DIR: gwas_hits.csv, causal_snps_truth.csv, tested.csv,
              detection_counts.csv, populations.csv (read by 06_analyze_gwas_results.R).
HOW TO RUN:   Rscript --vanilla R/scenarios/combine_scenarios.R   (from the repo root)
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  SCENARIO_ROOT  path, default output/populations/       folder with the scenario_* populations
  SCENARIOS      character, default all scenario_<h2>_<Mc> folders in SCENARIO_ROOT
  P_THRESHOLD    numeric, default 1e-3                   P value below which a SNP is a hit [COMMON]
  N_POINTS       integer, default 6                      resolution n_<N_POINTS> whose GWAS is used
  OUT_DIR        path, default output/combined_scenarios where the tables are written
=============================================================================
```

### R/analysis/regression_sumstats_dataprep.R

```text
=============================================================================
regression_sumstats_dataprep.R — PPIs of one population at all resolutions
=============================================================================
WHAT IT DOES: For every resolution n in RESOLUTIONS measures the true reaction norms of
  one simulated population at n equidistant points (same genotypes and plants; only
  the sampling changes) and computes all plasticity scores and the nine summary
  statistics. Every resolution is written as its own folder, exactly as when done one
  by one (sample_resolution(), compute_population_ppis()). A resolution that already
  has its PPIs is kept. These folders are the input of resolution_analysis.R,
  combined_resolution_figure.R, correlation_resolution.R and
  correlation_resolution_experiment.R.
REQUIRES:     POPULATION_DIR: a simulated population folder (simulate_population()).
PRODUCES:     POPULATION_DIR/n_<n>/ for every n: data/ (measurements), variance/,
              ppis/ (scores), analyses/summary_statistics.csv
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R","analysis","regression_sumstats_dataprep.R"))
-----------------------------------------------------------------------------
PARAMETERS — POPULATION_DIR is set before sourcing; the resolutions are edited here
  POPULATION_DIR  path, REQUIRED   simulated population folder                         [COMMON]
  RESOLUTIONS     integer, default c(2, 3, 4, 5, 6, 11, 26, 50)  sampling points        [COMMON]
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
REQUIRES:     RESOLUTION_DIR with its PPIs (compute_population_ppis()); the scores and
              summary statistics are those of the population's own sampling resolution.
PRODUCES:     Coefficient tables and diagnostic PDFs under RESOLUTION_DIR/plots/regression_summary_stats/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","analysis","regression_summarystats.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  RESOLUTION_DIR   path, REQUIRED (set before sourcing)   resolution folder               [COMMON]
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
WHAT IT DOES: Ranks the genotypes under every plasticity index and computes the
  Kendall rank correlation of every index pair, for the population at its own
  sampling resolution.
REQUIRES:     RESOLUTION_DIR with its PPIs (compute_population_ppis()).
PRODUCES:     Rank-correlation tables under RESOLUTION_DIR/analyses/rank_flops/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","analysis","kendall_rank_corr.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  RESOLUTION_DIR        path, REQUIRED                  resolution folder               [COMMON]
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
REQUIRES:     kendall_rank_corr.R run on RESOLUTION_DIR (analyses/rank_flops/); the
              summary statistics are those of the population at its own resolution.
PRODUCES:     Clustering summaries/plots of index agreement in
              RESOLUTION_DIR/analyses/rank_flops/derived/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","analysis","kendall_tau_corr_clustering.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  RESOLUTION_DIR  path, REQUIRED (set before sourcing)   resolution folder              [COMMON]
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
REQUIRES:     03 (so RESOLUTION_DIR must be set: a resolution folder of a simulated population).
PRODUCES:     Clustering figures under RESOLUTION_DIR/analyses/clustering/.
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
REQUIRES:     03 (so RESOLUTION_DIR must be set: a resolution folder of a simulated population).
PRODUCES:     Summary-stat vs score figures under RESOLUTION_DIR/analyses/summary_stats_vs_scores/.
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
WHAT IT DOES: Draws the sampling-sensitivity figure of a population: left, the slope
  of the regression of the PPI-summary-statistic correlation on the number of
  sampling points (heatmap, from R/analysis/resolution_analysis.R); right, the effect
  of the number of sampling points on every PPI from the mixed model
  PPI ~ n + summary statistics + (1 | genotype), fitted here on the resolution tables
  (forest plot).
REQUIRES:     The resolution folders POPULATION_DIR/n_* (regression_sumstats_dataprep.R)
              and POPULATION_DIR/analyses/resolution/ (resolution_analysis.R).
PRODUCES:     POPULATION_DIR/plots/combined_resolution_figure.pdf and the mixed-model
              results POPULATION_DIR/analyses/resolution/mixed_model.csv
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R","figures","combined_resolution_figure.R"))
-----------------------------------------------------------------------------
PARAMETERS — set before sourcing
  POPULATION_DIR  path, REQUIRED   population folder with the resolution analysis      [COMMON]
=============================================================================
```

### R/figures/correlation_resolution.R

```text
=============================================================================
correlation_resolution.R — correlation vs sampling resolution
=============================================================================
WHAT IT DOES: Plots how index correlations change with the sampling resolution of the
  reaction norm.
REQUIRES:     The resolution folders POPULATION_DIR/n_2 ... n_50 with their PPIs
              (regression_sumstats_dataprep.R).
PRODUCES:     boxplot_resolution.pdf, slope_heatmap_resolution.pdf in POPULATION_DIR/plots/.
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R","figures","correlation_resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder                [COMMON]
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
REQUIRES:     The resolution folders POPULATION_DIR/n_2 ... n_50 with their PPIs
              (regression_sumstats_dataprep.R).
PRODUCES:     confoundermodel_sumstats.pdf in POPULATION_DIR/plots/.
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R","figures","correlation_resolution_experiment.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder                [COMMON]
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
REQUIRES:     RESOLUTION_DIR with its PPIs (compute_population_ppis()); scores and
              summary statistics at the population's own sampling resolution.
PRODUCES:     Correlation-panel figures in RESOLUTION_DIR/plots/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","figures","correlation_summary_stats_PS.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  RESOLUTION_DIR  path, REQUIRED (set before sourcing)   resolution folder                [COMMON]
  stat_cols   character vector  summary-stat columns to use
=============================================================================
```

### R/figures/plasticity_score_realisation_plot.R

```text
=============================================================================
plasticity_score_realisation_plot.R — plasticity-score realisations (simulated)
=============================================================================
WHAT IT DOES: Plots the realised plasticity scores of a population at its own
  sampling resolution: distributions, Kendall-distance clustering per shape and
  combined, silhouette choice of k, Mantel tests between the shapes (Figure 2).
REQUIRES:     RESOLUTION_DIR with its PPIs (RESOLUTION_DIR/ppis/ppi_scores.rds from
              compute_population_ppis(RESOLUTION_DIR)).
PRODUCES:     Figures under RESOLUTION_DIR/plots/ (figure2_2.pdf, figure2_3.pdf, ...).
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","figures","plasticity_score_realisation_plot.R"))
-----------------------------------------------------------------------------
PARAMETERS — set before sourcing
  RESOLUTION_DIR  path, REQUIRED   resolution folder                                     [COMMON]
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
plotting_reaction_norms_all.R — reaction-norm overview figure of a population
=============================================================================
WHAT IT DOES: Draws the reaction norms of a population folder: per shape (or trait)
  selected genotypes with their replicate measurements (dashed) and replicate range
  (ribbon), plus histograms of the measured trait values. For a simulated population
  the solid line is the mother reaction norm (smooth curve from the mother
  parameters); for own data it is the genotype mean over the replicates.
REQUIRES:     A resolution folder (data/), simulated with simulate_population() /
              a scenario driver, or own data placed in data/.
PRODUCES:     RESOLUTION_DIR/plots/reaction_norms_final_figure.pdf,
              reaction_norms_final_figure_panels.pdf, trait_histogram_combined.pdf.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "my_population", "n_6")
              source(here::here("R","figures","plotting_reaction_norms_all.R"))
-----------------------------------------------------------------------------
PARAMETERS — set these before sourcing (guarded defaults)
  RESOLUTION_DIR  path, REQUIRED   resolution folder to draw                            [COMMON]
  PLOT_GENOTYPES  list shape/trait -> genotype IDs, default: 3 genotypes per panel at the
                  10/50/90 % quantile of the amplitude a (simulated) or of the genotype
                  mean trait value (own data)
  N_CURVE_POINTS  integer, default 200   points of the smooth mother curves (simulated)
=============================================================================
```

### R/figures/plotting_resolution.R

```text
=============================================================================
plotting_resolution.R — resolution correlation plots
=============================================================================
WHAT IT DOES: Reads the per-resolution correlation CSVs and renders the resolution
  correlation plots across the defined reaction-norm ranges.
REQUIRES:     correlations_*.csv of resolution.R in POPULATION_DIR/analyses/correlation_summary_stats/.
PRODUCES:     plots_<range>_<analysis>_per_score.pdf in POPULATION_DIR/plots/.
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R","figures","plotting_resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS (edit the named assignment near the top of the script)
  POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder                [COMMON]
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
  Includes the top-k overlap of the most plastic genotypes between indices (Figure 3).
REQUIRES:     Simulated data: RESOLUTION_DIR with its PPIs (ppis/ppi_scores.rds from
              compute_population_ppis(RESOLUTION_DIR)); the PPIs and summary statistics
              are those of the population's own sampling resolution.
              Maize (form_ranges of length 3): data/maize_scores_list.rds.
PRODUCES:     RESOLUTION_DIR/plots/plasticity_comprehensive_analysis/ (PDF and CSV tables);
              maize: output/plasticity_comprehensive_analysis_maize/.
HOW TO RUN:   setwd("~/PlastQuant")
              RESOLUTION_DIR <- here::here("output", "populations", "ppi_comparison", "n_6")
              source(here::here("R","figures","rank_corr_biplot.R"))
-----------------------------------------------------------------------------
PARAMETERS
  RESOLUTION_DIR  path, REQUIRED for simulated data                                       [COMMON]
=============================================================================
```

### R/figures/resolution.R

```text
=============================================================================
resolution.R — resolution experiment (dataprep + regression combined)
=============================================================================
WHAT IT DOES: Runs, in one pass, what regression_sumstats_dataprep.R and
  regression_summarystats.R do together, for the resolution experiment: Kendall
  correlations (permutation test) between every score and every summary statistic,
  pooled over the forms and per form, for every sub-range and sampling interval.
REQUIRES:     POPULATION_DIR: a simulated population folder; its reaction norms are
              sampled at every range x interval of a grid of global_initial_length points.
PRODUCES:     correlations_<range>_interval_<i>_indices_<...>.csv in
              POPULATION_DIR/analyses/correlation_summary_stats/ (read by plotting_resolution.R).
HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","resolution.R"))
-----------------------------------------------------------------------------
PARAMETERS: inherits 03's parameters via the auto-source; resolution knobs are set in
  the body of the script (edit them there).
=============================================================================
```

### R/analysis/resolution_analysis.R

```text
=============================================================================
resolution_analysis.R — sampling sensitivity: correlation and its regression
=============================================================================
WHAT IT DOES: From the resolution folders n_2 ... n_50 of a population (the same
  reaction norms measured at 2, 3, 4, 5, 6, 11, 26 and 50 points; their PPIs and
  summary statistics from regression_sumstats_dataprep.R) computes
    - the Pearson correlation of every PPI with every summary statistic across all
      genotypes, per number of sampling points,
    - for every PPI-statistic pair the linear regression of that correlation on the
      number of sampling points (slope, its standard error, t-test).
  The figure script R/figures/combined_resolution_figure.R reads these results and
  fits the mixed model on the tables.
REQUIRES:     regression_sumstats_dataprep.R run on POPULATION_DIR.
              The results span all resolutions and are written at population level.
PRODUCES:     POPULATION_DIR/analyses/resolution/
                correlations.csv            N_points, PPI, Stat, r
                correlation_regression.csv  PPI, Stat, intercept, slope, se_slope, t_stat, p_value, nobs
HOW TO RUN:   setwd("~/PlastQuant")
              POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
              source(here::here("R", "analysis", "resolution_analysis.R"))
-----------------------------------------------------------------------------
PARAMETERS
  POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder               [COMMON]
=============================================================================
```
