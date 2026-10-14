# PlastQuant

A reproducible simulation, GWAS, and figure pipeline for studying phenotypic
plasticity indices. PlastQuant simulates genotypes and reaction norms, computes a
broad panel of plasticity indices via the companion **`ppindices`** package, runs
GWAS on the resulting scores, and produces the manuscript figures.

Plasticity-index definitions live in the standalone **`ppindices`** package; this
repository is the *study* layer that depends on it. All paths are project-relative
(via [`here`](https://here.r-lib.org/)) and the R environment is locked with
[`renv`](https://rstudio.github.io/renv/).

## Requirements

- R (>= 4.1)
- A C/C++ toolchain (some dependencies, e.g. `sommer`, compile from source)
- The **`ppindices`** companion package. It is pinned in `renv.lock` and installed
  by `renv::restore()`; make sure its source is available locally (as a sibling
  `~/ppindices`) before restoring.

## Setup

```r
setwd("~/PlastQuant")        # always work from the project root
install.packages("renv")     # if not already available
renv::restore()              # installs the exact package versions from renv.lock
```

`renv::restore()` provisions everything, including `ppindices`, into a project-local
library. Nothing is installed into your system library.

> **Tip:** Start R from the project root (or open the project) so `here::here()`
> resolves to the repository, and use a **clean session** — do not restore an old
> workspace (`R --no-restore-data`, or turn off "Restore .RData" in RStudio). A
> restored workspace can carry stale variables from a previous run and misdirect
> outputs.

## Repository layout

```
R/                                pipeline, analysis, figure, and scenario scripts
  sim/                            simulation model: SNPs -> genetic values -> reaction norms
                                  -> measurements; all variance settings in sim_config()
  ppi/                            PPI registry and computation (uses ppindices)
  pipeline/                       input contract, GWAS, power per shape, pipeline runner
  03_plasticity_scores.R          simulate (R/sim) and compute all plasticity indices
  03_plasticity_scores_maize.R    same, for the empirical maize dataset
  04_run_gwas.R                   GWAS on the plasticity scores
  05_evaluate_gwas.R              power per shape of one scenario run (heatmaps)
  06_analyze_gwas_results.R       summarise the GWAS results of all scenarios
  06_analyze_maize_gwas_results.R summarise maize-GWAS results
  scenarios/                      per-scenario driver scripts, scenario_maize.R,
                                  combine_scenarios.R (collects all scenario runs)
  analysis/                       regression / rank-correlation / clustering summaries
  figures/                        figure-producing scripts
data/                             committed input datasets
  A-Datasets/                     empirical / maize input data files
  maize_scores_list.rds           precomputed maize plasticity scores
output/                           script outputs
  regression_summary_stats/       committed — regression inputs for figures
  correlation_summary_stats/      committed — correlation inputs for figures
  rank_flops/                     committed — rank-agreement inputs for figures
  scenario_*/  plots/  ...        generated at run time (git-ignored)
tests/                            unit tests (Rscript --vanilla tests/run_tests.R)
docs/PIPELINE.md                  detailed per-script reference
```

Every script begins with a header block documenting what it does, what it requires,
what it produces, how to run it, and its parameters. **[`docs/PIPELINE.md`](docs/PIPELINE.md)**
collects all of these plus a cross-script dependency map — start there for detail.

## Running the pipeline

Each **scenario** driver is a single call that simulates the population, computes all
PPIs, runs the GWAS (rrBLUP mixed model with kinship) and evaluates the power per
shape. From the repository root:

```bash
Rscript --vanilla R/scenarios/scenario_20_5.R
```

Scenario drivers are named `scenario_<h2>_<causal_snps>.R`. The first token is h² in
percent (`20`/`40`/`60` → `h2 = 0.2/0.4/0.6`); the second is the number of causal SNPs
per parameter. Each driver calls

```r
source(here::here("R", "pipeline", "load_pipeline.R"))
run_simulation_scenario(sim_config(h2 = 0.2, Mc = 5), out_dir = here::here("output", "scenario_20_5"))
```

All other settings and their defaults are arguments of `sim_config()` in
`R/sim/simulation_config.R`:

- variance components (`sigma_L`, `repeatability_beta`, `error_share`)
- number of genotypes, replicates, SNPs and measurement points
- seed and study space

The same seed gives the same genetics for every h². Outputs are written to
`output/scenario_<...>/run_<timestamp>/`, including the realised variance components
(`variance_latent.csv`, `variance_trait.csv`). Each driver then sources
`05_evaluate_gwas.R`, which writes the power per shape and parameter
(`gwas_accuracy_summary.csv`) and the heatmaps (`plots/gwas_power_heatmap.pdf`). A
parameter's power is always computed within its shape. A hit on a causal SNP of another
shape counts as a false positive.

To compare the scenarios, collect their newest runs and summarise them:

```bash
Rscript --vanilla R/scenarios/combine_scenarios.R    # tables in output/combined_scenarios/
Rscript --vanilla R/06_analyze_gwas_results.R        # plots and precision/recall tables
```

`scenario_maize.R` runs the GWAS on the empirical maize scores.

**Own data.** Pass a SNP matrix, reaction norms and, optionally, known causal SNPs:

```r
source(here::here("R", "pipeline", "load_pipeline.R"))
input <- pipeline_input(phenotypes, environments, snp_matrix = X, truth = known_snps)
run_pipeline(input, out_dir = here::here("output", "my_data"))
```

Without SNPs, use `run_pipeline(pipeline_input(phenotypes, environments), out_dir, gwas = FALSE)`;
this computes the PPIs only.

**PPIs.** To add a PPI, register it:

```r
reg <- register_ppi(default_ppi_registry(), ppi_entry("Name", fun, inputs = c(trait_values = "trait", env_values = "env")))
```

Pass `registry = reg` to the runner. To leave PPIs out, use `exclude = c(...)`. FW is
the one PPI computed on the whole population's reaction norms. It is always computed
once on the replicate means and gives one value per genotype.

See `docs/PIPELINE.md` for the data format, the output files and the numbered scripts.

## Reproducing the figures

The regression, correlation, and rank-agreement datasets the figures consume are
**committed under `output/`**, so most figure scripts run directly, without re-running
the whole pipeline. Figure PDFs are written to `output/plots/`.

```r
setwd("~/PlastQuant")
source(here::here("R", "figures", "plasticity_score_realisation_plot.R"))
```

| Figure / artefact                         | Script (`R/figures/`)                       | Needs a scenario run first? |
|-------------------------------------------|---------------------------------------------|-----------------------------|
| Reaction-norm overview                    | `plotting_reaction_norms_all.R`             | yes — run `scenario_20_5.R` |
| Plasticity-score realisations (simulated) | `plasticity_score_realisation_plot.R`       | no                          |
| Plasticity-score realisations (maize)     | `plasticity_score_realisation_plot_maize.R` | no                          |
| Summary-stats correlation panels          | `correlation_summary_stats_PS.R`            | no                          |
| Resolution / correlation experiments      | `correlation_resolution.R`, `correlation_resolution_experiment.R`, `combined_resolution_figure.R`, `plotting_resolution.R` | no |
| Rank-correlation biplot                   | `rank_corr_biplot.R`                         | no                          |

The committed datasets can be regenerated: `R/analysis/regression_sumstats_dataprep.R`
rebuilds `output/regression_summary_stats/`, and the Kendall / resolution scripts
(`R/analysis/kendall_*`, `R/figures/resolution.R`) rebuild the correlation and
rank-agreement CSVs.

## Notes

- **Warnings are suppressed** project-wide (`options(warn = -1)` in `.Rprofile` and at
  the top of every script) to keep console output readable. Set `options(warn = 0)` in
  your session if you need to see them while debugging.

## Dependencies

The full, version-pinned dependency set is recorded in `renv.lock`. Key packages:
`ppindices`, `here`, `rrBLUP`, `sommer` (GWAS); `dplyr`, `tidyr`, `purrr`, `readr`,
`reshape2` (data wrangling); `ggplot2`, `gridExtra`, `patchwork`, `pheatmap`,
`ggpubr` (figures); `vegan`, `dendextend`, `factoextra`, `mclust` (multivariate /
clustering analyses).

## License

MIT.
