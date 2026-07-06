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
  01_simulate_genetics.R          simulate SNP genotypes + genetic parameters
  02_reaction_norms_linear.R      generate linear reaction norms
  02_reaction_norms_nonlinear.R   generate Gaussian/sinusoidal/wave reaction norms
  03_plasticity_scores.R          compute all plasticity indices (uses ppindices)
  03_plasticity_scores_maize.R    same, for the empirical maize dataset
  04_run_gwas.R                   GWAS on the plasticity scores
  05_evaluate_gwas.R              accuracy evaluation of GWAS hits
  06_analyze_gwas_results.R       summarise simulated-GWAS results
  06_analyze_maize_gwas_results.R summarise maize-GWAS results
  scenarios/                      per-scenario driver scripts + scenario_maize.R
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
docs/PIPELINE.md                  detailed per-script reference
```

Every script begins with a header block documenting what it does, what it requires,
what it produces, how to run it, and its parameters. **[`docs/PIPELINE.md`](docs/PIPELINE.md)**
collects all of these plus a cross-script dependency map — start there for detail.

## Running the pipeline

Each **scenario** driver sets its parameters and runs the full
simulate → scores → GWAS → evaluate chain. From the repository root:

```r
setwd("~/PlastQuant")
source(here::here("R", "scenarios", "scenario_20_5.R"))
```

Scenario drivers are named `scenario_<H>_<causal_snps>.R`, where the first token maps
to the `HERITABILITY_5TH` setting (`20`/`40`/`60` → `1`/`2`/`3`) and the second is the
number of causal SNPs per parameter. `scenario_maize.R` runs the empirical maize
analysis. Outputs are written under `output/scenario_<...>/`.

To change a parameter, edit the clearly-marked block at the top of the relevant script
(or set it in your session before sourcing). You never edit the pipeline internals —
see each script's header, or `docs/PIPELINE.md`, for the full list of parameters.

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

> Note: `R/scenarios/combine_scenarios_snps.sh` aggregates per-scenario GWAS CSVs. It
> uses relative paths and expects to be run from inside `output/`, where the
> `scenario_*` folders live.

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
