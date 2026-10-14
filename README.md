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
>
> If `here::here()` points somewhere else (e.g. `cannot open file '~/R/pipeline/...'`),
> R was started outside the project. Restart R in `~/PlastQuant`, or run
> `setwd("~/PlastQuant"); here::i_am("README.md")`.

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

The pipeline has four steps. Each step is a function and can be run on its own:

    simulation ─▶ PPIs ─▶ GWAS ─▶ evaluation (power per shape) ─▶ plots
    simulate_reaction_norms()  compute_ppis()  run_gwas()  evaluate_power()

Load the functions once per R session, from the project root:

```r
source(here::here("R", "pipeline", "load_pipeline.R"))
```

### 1. Simulation only

```r
cfg <- sim_config(h2 = 0.4, Mc = 5)   # every setting is described in R/sim/simulation_config.R
sim <- simulate_reaction_norms(cfg)
```

`sim_config()` takes these settings:

- `h2`: share of the latent genetic variance from the causal SNPs.
- `Mc`: causal SNPs per parameter.
- `n_per_shape` (200), `n_rep` (3), `n_snps` (200000), `n_points` (6): population size and design.
- `seed`.
- `sigma_L`, `repeatability_beta`, `error_share`: the variance components.
- `study_space`: parameter ranges.

The same `seed` gives the same genetics for every `h2`.

`sim` contains:

| Element | Content |
|---|---|
| `sim$phenotypes` | list with one table per shape (`linear`, `gaussian`, `sinusoidal`, `wave`); columns `Genotype`, `Replicate`, `Environment`, `Trait` (the measured values) |
| `sim$environments` | environment values of the measurement points (`1 + 9 t`) |
| `sim$snp_matrix` | genotypes × SNPs, dosages 0/1/2, row names = genotype IDs |
| `sim$truth` | causal SNPs: `Shape`, `Parameter`, `SNP`, `Effect` |
| `sim$mother_parameters` | parameter values of every genotype's mother reaction norm |
| `sim$mother_curves` | mother reaction norms at the measurement points |
| `sim$plant_parameters` | parameter values of every replicate plant |

To see the realised variance components:

```r
vr <- variance_report(sim)
vr$latent   # SNP share per parameter (equals h2) and mean repeatability
vr$trait    # trait level: causal SNPs, background, interaction, within genotype, measurement error
```

To keep the simulation, use `saveRDS(sim, here::here("output", "my_sim.rds"))`.

### 2. PPIs

```r
traits <- compute_ppis(sim$phenotypes, sim$environments)                         # all 28 PPIs
traits <- compute_ppis(sim$phenotypes, sim$environments, exclude = c("FW", "PPF")) # leave some out
traits$RN   # one table per PPI: Score, Type (shape), Genotype_ID, Replicate
```

By default, the replicates are averaged per genotype and environment before a PPI is
computed (`replicate_handling = "mean_before"`). With `replicate_handling = "replicates"`,
each replicate reaction norm gets its own PPI value (one row per replicate). FW is the
exception: it uses the reaction norms of the whole population, so it is always
computed once on the replicate means, with one value per genotype.

To add your own PPI:

```r
reg <- register_ppi(default_ppi_registry(),
                    ppi_entry("MyIndex", function(trait_values, env_values) diff(range(trait_values)) / diff(range(env_values)),
                              inputs = c(trait_values = "trait", env_values = "env")))
traits <- compute_ppis(sim$phenotypes, sim$environments, registry = reg)
```

### 3. GWAS and power

```r
gwas  <- run_gwas(traits, sim$snp_matrix)          # rrBLUP mixed model with kinship, per PPI and shape
power <- evaluate_power(gwas, sim$truth, n_snps = ncol(sim$snp_matrix))   # Bonferroni 0.05 / n_SNPs
```

`power` has one row per shape, PPI and parameter. A parameter is only evaluated within its
own shape, and a hit on a SNP that is not causal for that shape counts as a false positive.

### Steps 2–3 in one call, with all files written

```r
run_dir <- run_pipeline(sim, out_dir = here::here("output", "my_simulation"))
```

For a simulation, use `run_simulation_scenario()` instead. It simulates and then calls
`run_pipeline()`. It also writes the variance report, the mother parameters and the
settings (`config.rds`):

```r
run_dir <- run_simulation_scenario(sim_config(h2 = 0.3, Mc = 8), out_dir = here::here("output", "my_simulation"))
```

A run folder `out_dir/run_<timestamp>/` contains:

- `ppi_scores/<PPI>.csv` and `ppi_scores.rds`
- `phenotypes.rds` and `snp_matrix.rds`
- `truth.csv`
- `gwas_results.csv` and `gwas_power.csv`
- for scenarios also `variance_latent.csv`, `variance_trait.csv`, `measurement.csv`,
  `mother_parameters.csv`, `mother_curves.rds` and `config.rds`

### 4. Plots

- **Reaction norms** of a simulation: mother curves, replicates and trait histograms.

  ```r
  sim <- simulate_reaction_norms(sim_config(h2 = 0.4, Mc = 5, n_snps = 2000))   # or readRDS(...)
  source(here::here("R", "figures", "plotting_reaction_norms_all.R"))           # draws `sim`
  ```

  Without a `sim` in the session, the script simulates its own population.

- **Power heatmaps** per shape for one run:

  ```r
  RUN_DIR <- run_dir
  source(here::here("R", "05_evaluate_gwas.R"))   # RUN_DIR/plots/gwas_power_heatmap.pdf
  ```

- **Comparison of the heritability scenarios**: see below.

### Heritability scenarios

The nine drivers in `R/scenarios/` simulate the heritability grid h² ∈ {0.2, 0.4, 0.6} ×
Mc ∈ {5, 10, 20}. Each driver is `run_simulation_scenario()` with its h² and Mc,
followed by 05:

```bash
Rscript --vanilla R/scenarios/scenario_20_5.R        # h2 = 0.2, 5 causal SNPs per parameter
```

After running the scenarios, collect their newest runs and summarise them:

```bash
Rscript --vanilla R/scenarios/combine_scenarios.R    # tables in output/combined_scenarios/
Rscript --vanilla R/06_analyze_gwas_results.R        # plots, precision/recall tables
```

### Own data

Your data goes into the same steps. You need:

- **Phenotypes**: a named list with one table per trait (or shape), columns `Genotype`,
  `Replicate`, `Environment`, `Trait`. Use `Replicate = 1` if there are no replicates.
- **Environments**: the numeric environment values, sorted. Every `Environment` in the
  tables must be one of them.
- **SNP matrix** (optional): genotypes × SNPs, dosages 0/1/2, with the genotype IDs of the
  phenotypes as row names and SNP IDs as column names.
- **Truth** (optional, needs SNPs): known causal SNPs, e.g. known plasticity genes.
  Columns `Shape` (= a name of the phenotype list), `Parameter` (free label) and `SNP`.

From a long CSV (one row per measurement):

```r
d <- read.csv(here::here("data", "my_phenotypes.csv"))     # Genotype, Replicate, Environment, Trait
phenotypes   <- list(biomass = d[, c("Genotype", "Replicate", "Environment", "Trait")])
environments <- sort(unique(d$Environment))
```

From a wide table (one column per environment):

```r
wide <- read.csv(here::here("data", "my_wide.csv"))       # Genotype, Replicate, E1, E5, E9, ...
envs <- c(E1 = 1, E5 = 5, E9 = 9)                          # column name -> environment value
long <- do.call(rbind, lapply(names(envs), function(e)
  data.frame(Genotype = wide$Genotype, Replicate = wide$Replicate, Environment = envs[[e]], Trait = wide[[e]])))
phenotypes   <- list(height = long)
environments <- unname(envs)
```

**PPIs only (no SNPs):**

```r
input   <- pipeline_input(phenotypes, environments)        # checks the format, names what is wrong
traits  <- compute_ppis(input$phenotypes, input$environments)
run_dir <- run_pipeline(input, out_dir = here::here("output", "my_data"), gwas = FALSE)   # writes ppi_scores/
```

**PPIs and GWAS:**

```r
snps  <- as.matrix(read.csv(here::here("data", "my_snps.csv"), row.names = 1, check.names = FALSE))
input <- pipeline_input(phenotypes, environments, snp_matrix = snps)
run_dir <- run_pipeline(input, out_dir = here::here("output", "my_data"))
```

**With known causal SNPs**, power is computed as well:

```r
truth <- data.frame(Shape = "biomass", Parameter = "known", SNP = c("SNP_12", "SNP_80"))
input <- pipeline_input(phenotypes, environments, snp_matrix = snps, truth = truth)
run_dir <- run_pipeline(input, out_dir = here::here("output", "my_data"))
RUN_DIR <- run_dir; source(here::here("R", "05_evaluate_gwas.R"))   # power heatmap
```

Notes on own data:

- Some PPIs need a minimum number of environments. RNN uses a degree-3 polynomial and
  needs more than 3. A PPI that cannot be computed gives `NA` with a message; the run
  continues.
- Genotype IDs can be text or numbers.

### What works with which data

| Plot / table | Simulated data | Own data |
|---|---|---|
| PPI tables (`ppi_scores/*.csv`) | yes | yes |
| GWAS results, power table | yes | GWAS: with SNPs; power: with truth |
| Power heatmaps (`05_evaluate_gwas.R`) | yes | with SNPs and truth |
| Reaction-norm figure (`plotting_reaction_norms_all.R`) | yes | no (needs mother parameters) |
| Scenario comparison (`combine_scenarios.R`, `06_analyze_gwas_results.R`) | scenario runs | no |
| PPI distribution / clustering figures (`plasticity_score_realisation_plot.R`, `rank_corr_biplot.R`) | not yet: they read the committed regression table or the maize scores | not yet |

### The scores script and the analyses

`R/03_plasticity_scores.R` simulates a population with the same model and computes all
PPIs into `scores_list`, the input of the analyses in `R/analysis/` and of
`R/figures/resolution.R`. Those scripts source 03 themselves. To run 03 directly:

```r
H2 <- 0.4; CAUSAL_SNP_NUM <- 5; OUTPUT_BASE <- here::here("output", "my_run")
source(here::here("R", "03_plasticity_scores.R"))
```

`scenario_maize.R` runs the GWAS on the empirical maize scores.

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
| Reaction-norm overview                    | `plotting_reaction_norms_all.R`             | no (simulates its own population) |
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
