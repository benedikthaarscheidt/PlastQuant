# PlastQuant

A reproducible simulation, GWAS, and figure pipeline for studying phenotypic
plasticity indices. PlastQuant simulates genotypes and reaction norms, computes a
broad panel of plasticity indices via the [`ppindices`](https://…) package, runs
GWAS on the resulting scores, and produces the manuscript figures.

Plasticity-index definitions live in the standalone **`ppindices`** package; this
repository is the *study* layer that depends on it. All paths are project-relative
(via [`here`](https://here.r-lib.org/)) and the R environment is locked with
[`renv`](https://rstudio.github.io/renv/).

## Requirements

- R (>= 4.1)
- The `ppindices` package (installed into the project library by `renv::restore()`)
- Pandoc (only if you rebuild R Markdown outputs)

## Setup

```r
# from the repository root
install.packages("renv")     # if not already available
renv::restore()              # installs the exact package versions from renv.lock
```

`renv::restore()` provisions everything, including `ppindices`, into a project-local
library. Nothing is installed into your system library.

## Repository layout

```
R/                      numbered pipeline (run in order)
  01_simulate_genetics.R          simulate SNP genotypes + genetic parameters
  02_reaction_norms_linear.R      generate linear reaction norms
  02_reaction_norms_nonlinear.R   generate Gaussian/sinusoidal/wave reaction norms
  03_plasticity_scores.R          compute all plasticity indices (uses ppindices)
  03_plasticity_scores_maize.R    same, for the empirical maize dataset
  04_run_gwas.R                   GWAS on the plasticity scores
  05_evaluate_gwas.R              accuracy evaluation of GWAS hits
  06_analyze_gwas_results.R       summarise simulated-GWAS results
  06_analyze_maize_gwas_results.R summarise maize-GWAS results
  scenarios/            per-scenario driver scripts (scenario_<geno>_<snp>.R)
  analysis/             regression / rank-correlation / clustering summaries
  figures/              figure-producing scripts
data/                   input datasets (A-Datasets/, maize_scores_list.rds)
output/                 all generated artefacts (git-ignored)
```

## Running the pipeline

Each **scenario** driver sets its parameters and sources the pipeline. From the
repository root:

```r
setwd("~/PlastQuant")   # or open the project so here::here() resolves to the root
source(here::here("R", "scenarios", "scenario_20_5.R"))
```

Scenario drivers are named `scenario_<num_genotypes>_<num_causal_snps>.R`
(e.g. `scenario_20_5.R` = smallest). `scenario_maize.R` runs the empirical maize
analysis. Outputs are written under `output/scenario_<...>/`.

## Reproducing the figures

Figure scripts read the pipeline outputs and write PDFs into
`output/plots/figures/`. Run the relevant scenarios first, then:

| Figure / artefact                         | Script (`R/figures/`)                       |
|-------------------------------------------|---------------------------------------------|
| Reaction-norm overview                    | `plotting_reaction_norms_all.R`             |
| Plasticity-score realisations (simulated) | `plasticity_score_realisation_plot.R`       |
| Plasticity-score realisations (maize)     | `plasticity_score_realisation_plot_maize.R` |
| Summary-stats correlation panels          | `correlation_summary_stats_PS.R`            |
| Resolution / correlation experiments      | `correlation_resolution.R`, `correlation_resolution_experiment.R`, `combined_resolution_figure.R` |
| Rank-correlation biplot                   | `rank_corr_biplot.R`                         |

Analysis summaries used as figure inputs are produced by the scripts in
`R/analysis/` (regression summary stats, Kendall rank correlations, clustering).

> Note: `R/scenarios/combine_scenarios_snps.sh` is a manual helper that aggregates
> per-scenario GWAS CSVs. It uses relative paths and expects to be run from the
> directory that holds the `scenario_*` output folders (i.e. inside `output/`).

## Dependencies

The full, version-pinned dependency set is recorded in `renv.lock`. Key packages:
`ppindices`, `here`, `rrBLUP`, `sommer` (GWAS); `dplyr`, `tidyr`, `purrr`, `readr`,
`reshape2` (data wrangling); `ggplot2`, `gridExtra`, `patchwork`, `pheatmap`,
`ggpubr` (figures); `vegan`, `dendextend`, `factoextra`, `mclust` (multivariate /
clustering analyses).

## License

MIT.
