# PlastQuant

A reproducible simulation, GWAS, and figure pipeline for studying phenotypic
plasticity indices. PlastQuant simulates genotypes and their reaction norms, measures
them at chosen sampling resolutions, computes a broad panel of plasticity indices
(PPIs) via the companion **`ppindices`** package, runs GWAS on the resulting scores,
and produces the manuscript figures.

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
> workspace (`R --no-restore-data`, or turn off "Restore .RData" in RStudio).
>
> If `here::here()` points somewhere else (e.g. `cannot open file '~/R/pipeline/...'`),
> R was started outside the project. Restart R in `~/PlastQuant`, or run
> `setwd("~/PlastQuant"); here::i_am("README.md")`.

## Repository layout

```
R/
  pipeline/                       the infrastructure: load_pipeline.R (loads all functions),
                                  population_folder.R (folder format), pipeline_input.R
                                  (data contract), sample_resolution.R (measuring at a resolution)
  sim/                            the simulation model (SNPs -> genetic values -> true reaction
                                  norms; all settings in sim_config()), simulate_population.R,
                                  plot of the true norms
  ppi/                            PPI registry and computation (uses ppindices),
                                  compute_population_ppis.R, 03_plasticity_scores(_maize).R
  gwas/                           GWAS functions, run_population_gwas.R, power per shape,
                                  04_run_gwas.R, 05_evaluate_gwas.R (power heatmaps),
                                  06_analyze_gwas_results.R, 06_analyze_maize_gwas_results.R
    scenarios/                    heritability scenarios, combine_scenarios.R, scenario_maize.R
  analysis/                       comparison of the PPIs: analyses and their figures,
                                  resolution stability, summary statistics, rank correlation
data/                             committed input datasets (maize data, maize scores)
output/                           script outputs (git-ignored)
  populations/<name>/             one folder per population (see below)
tests/                            unit tests (Rscript --vanilla tests/run_tests.R)
docs/PIPELINE.md                  detailed per-script reference
```

Every script begins with a header block documenting what it does, what it requires,
what it produces, how to run it, and its parameters. **[`docs/PIPELINE.md`](docs/PIPELINE.md)**
collects all of these.

## The population folder

Everything belongs to a population and is written into its folder. The population
itself is stored once: for a simulation the true mother and replicate reaction norms,
without any measurement. Every **resolution** — the population measured at n
environments — is its own subfolder `n_<n>`, with the measurements and everything
computed from them.

```
output/populations/<name>/
  data/                    the population
    norms.rds                true mother and replicate reaction norms (all parameters)
    mother_parameters.csv    parameters of every genotype's mother norm (and repeatability)
    plant_parameters.csv     parameters of every replicate plant
    reaction_norms.pdf       plot of the true norms
    snp_matrix.rds|.csv      SNP dosages 0/1/2
    kinship.csv              optional; otherwise computed from the SNPs
    truth.csv                causal SNPs: Trait, Parameter, SNP
    config.rds               all simulation settings
  variance/                latent variance components (SNP share = h2, repeatability)
  n_6/                     the population measured at 6 points
    data/traits/<shape>.csv  measured values: Genotype, then one column per environment
    data/environments.csv    numeric value of every environment column
    variance/                trait-level variance shares and measurement error
    ppis/                    plasticity indices
    gwas/                    GWAS results and power
    analyses/  plots/        analyses and figures of this resolution
  n_11/ ...                same population, other resolutions
  analyses/  plots/        analyses and figures across resolutions
```

## Running the pipeline

Load the functions once per R session, from the project root (again after the code
has changed):

```r
source(here::here("R", "pipeline", "load_pipeline.R"))
```

### 1. Simulate a population

```r
POPULATION_DIR <- here::here("output", "populations", "My_population")
simulate_population(sim_config(h2 = 0.4, Mc = 5, n_per_shape = 20), POPULATION_DIR)
```

This simulates the SNPs, the genetic values and the true reaction norms of every
genotype (mother norm) and its replicate plants, for the four shapes (linear, Gaussian,
sinusoidal, wave), and writes `data/` (including the plot `data/reaction_norms.pdf`)
and `variance/`. Every setting of `sim_config()` is commented in
`R/sim/simulation_config.R`:

- `h2`: share of the latent genetic variance from the causal SNPs
- `Mc`: causal SNPs per parameter
- `n_per_shape` (200), `n_rep` (3), `n_snps` (200000): population size
- `seed` (1): the same seed gives the same genetics for every `h2`
- `sigma_L`, `repeatability_beta`, `error_share`: variance components
- `study_space`: parameter ranges

### 2. Measure it at a resolution

```r
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # -> n_6/
```

The true norms are measured at 6 equidistant environments (1 to 10), and the
measurement error is added (5 % of the trait variance at these points). The same plants
are measured at every resolution; measuring the same resolution again gives exactly the
same values.

### 3. PPIs

```r
compute_population_ppis(RESOLUTION_DIR)                       # all 28 PPIs -> n_6/ppis/
compute_population_ppis(RESOLUTION_DIR, exclude = c("PPF"))   # leave PPIs out
```

The argument `replicate_handling` decides on which reaction norms the PPIs are computed:

```r
compute_population_ppis(RESOLUTION_DIR, replicate_handling = "mean_before")  # default
compute_population_ppis(RESOLUTION_DIR, replicate_handling = "replicates")
```

- `"mean_before"` (the default, used when the argument is left out): the replicates of a
  genotype are averaged per environment first, and the PPI is computed on this averaged
  reaction norm — one value per genotype.
- `"replicates"`: the PPI is computed on every replicate's own measured reaction norm,
  nothing is averaged — one value per replicate (in the GWAS, the replicates of a
  genotype are then repeated observations of the same line).

FW uses the reaction norms of the whole population and is always computed once on the
averaged norms, one value per genotype, with either setting.

To add a PPI, register it and pass the registry:

```r
reg <- register_ppi(default_ppi_registry(),
                    ppi_entry("MyIndex", function(trait_values, env_values) diff(range(trait_values)) / diff(range(env_values)),
                              inputs = c(trait_values = "trait", env_values = "env")))
compute_population_ppis(RESOLUTION_DIR, registry = reg)
```

### 4. GWAS (optional)

```r
run_population_gwas(RESOLUTION_DIR)                  # -> n_6/gwas/
source(here::here("R", "gwas", "05_evaluate_gwas.R"))        # power heatmaps -> n_6/plots/
```

rrBLUP mixed model with the kinship matrix, one GWAS per PPI and shape. With a truth
table the power is computed per shape: a parameter is evaluated within its own shape,
and a hit on a SNP that is not causal for that shape counts as a false positive.

### Figures and analyses of one resolution

Set `RESOLUTION_DIR` (step 2) and source the script; it writes into that resolution's
`plots/` or `analyses/`. Steps 1–3 are needed.

```r
source(here::here("R", "analysis", "plasticity_score_realisation_plot.R"))  # Fig. 2: distributions, clustering, Mantel
source(here::here("R", "analysis", "rank_corr_biplot.R"))                   # Fig. 3: rank agreement, top-k overlap
source(here::here("R", "analysis", "plotting_reaction_norms_all.R"))        # measured reaction norms and histograms
source(here::here("R", "analysis", "correlation_summary_stats_PS.R"))       # PPIs vs summary statistics
source(here::here("R", "analysis", "regression_summarystats.R"))           # PPIs regressed on summary statistics
source(here::here("R", "analysis", "kendall_rank_corr.R"))                 # rank tables ...
source(here::here("R", "analysis", "kendall_tau_corr_clustering.R"))       # ... and their agreement analyses
source(here::here("R", "analysis", "clustering_rn_summary_stats.R"))       # clustering (via 03)
```

### Stability across resolutions

These scripts take `POPULATION_DIR`. The first creates the resolutions
n = 2, 3, 4, 5, 6, 11, 26, 50 of the population — steps 2 and 3 for each, exactly as by
hand — and saves the summary statistics of each; it is the slow part. A resolution that
already has its PPIs is kept.

```r
source(here::here("R", "analysis", "regression_sumstats_dataprep.R"))  # n_2 ... n_50
source(here::here("R", "analysis", "resolution_analysis.R"))           # correlation PPI-statistic per n and its
                                                                         # regression on n -> analyses/resolution/
source(here::here("R", "analysis", "combined_resolution_figure.R"))     # Fig. 4 (incl. mixed model) -> plots/
source(here::here("R", "analysis", "correlation_resolution.R"))
source(here::here("R", "analysis", "correlation_resolution_experiment.R"))
source(here::here("R", "analysis", "resolution.R"))                     # correlations on the gradient halves ...
source(here::here("R", "analysis", "plotting_resolution.R"))            # ... and their plots
```

### Heritability scenarios (GWAS power)

The nine drivers in `R/gwas/scenarios/` simulate the grid h² ∈ {0.2, 0.4, 0.6} ×
Mc ∈ {5, 10, 20} and run steps 1–4 at 6 points each:

```bash
Rscript --vanilla R/gwas/scenarios/scenario_20_5.R        # -> output/populations/scenario_20_5/
```

After running the scenarios, compare them:

```bash
Rscript --vanilla R/gwas/scenarios/combine_scenarios.R    # -> output/combined_scenarios/
Rscript --vanilla R/gwas/06_analyze_gwas_results.R        # plots, precision/recall tables
```

`scenario_maize.R` runs the GWAS on the empirical maize scores.

## Own data

Put your data into a population folder; the measurements form one resolution
`n_<number of environments>`:

```
output/populations/my_data/
  data/
    snp_matrix.csv      Genotype, SNP_1, SNP_2, ...   dosages 0/1/2        (optional)
    kinship.csv         Genotype x Genotype                                 (optional)
    truth.csv           Trait, Parameter, SNP: known causal SNPs            (optional)
  n_5/
    data/
      traits/biomass.csv  Genotype, E1, E2, E3, E4, E5   one row per plant: a genotype
                                                          appears once per replicate
      environments.csv    Column, Value                  numeric value of every column
```

Then run steps 3 and 4 and the figures of one resolution on it:

```r
RESOLUTION_DIR <- here::here("output", "populations", "my_data", "n_5")
compute_population_ppis(RESOLUTION_DIR)
run_population_gwas(RESOLUTION_DIR)        # needs snp_matrix.csv
```

Notes:

- Own data is measured at its environments only, so it cannot be measured at other
  resolutions; the stability analysis across resolutions needs a simulated population.
- Without `truth.csv` the GWAS gives results but no power.
- Some PPIs need a minimum number of environments (RNN: more than 3). A PPI that cannot
  be computed gives `NA` with a message; the run continues.
- Format errors are reported with the file and the offending column or genotype.

| Script | Simulated | Own data |
|---|---|---|
| PPIs, GWAS (`compute_population_ppis`, `run_population_gwas`) | yes | yes (GWAS with SNPs) |
| Power heatmaps (`05_evaluate_gwas.R`) | yes | with SNPs and truth |
| Fig. 2, Fig. 3, summary-statistic and Kendall analyses | yes | yes |
| Measured reaction norms (`plotting_reaction_norms_all.R`) | yes | yes (genotype means instead of mother norms) |
| `03_plasticity_scores.R` and the analyses that source it | yes | no (needs the four shapes) |
| Stability across resolutions | yes | no |
| Heritability scenarios | yes | no |

## Working in memory

The steps are thin wrappers around functions that can also be used directly:

```r
cfg      <- sim_config(h2 = 0.4, Mc = 5)
norms    <- simulate_norms(cfg)                         # true norms, SNPs, truth
measured <- measure_norms(norms, 6)                     # measurement at 6 points
traits   <- compute_ppis(measured$phenotypes, measured$environments)
gwas     <- run_gwas(traits, measured$snp_matrix)
power    <- evaluate_power(gwas, measured$truth, n_snps = ncol(measured$snp_matrix))
vr       <- variance_report(measured)                   # latent and trait-level variance
```

## Earlier results

`output/regression_summary_stats/`, `output/correlation_summary_stats/` and
`output/rank_flops/` hold committed tables of the earlier simulation. No script writes
or reads them any more; all scripts work on population folders.

## Notes

- **Warnings are suppressed** project-wide (`options(warn = -1)` in `.Rprofile` and at
  the top of every script) to keep console output readable. Set `options(warn = 0)` in
  your session if you need to see them while debugging.
- `README.md` is plain Markdown; GitHub and editors show it directly. Rendering it with
  `rmarkdown` additionally needs the program Pandoc (bundled with RStudio).

## Dependencies

The full, version-pinned dependency set is recorded in `renv.lock`. Key packages:
`ppindices`, `here`, `rrBLUP`, `sommer`, `lme4`, `lmerTest` (models); `dplyr`, `tidyr`,
`purrr`, `readr`, `reshape2` (data wrangling); `ggplot2`, `gridExtra`, `patchwork`,
`pheatmap`, `ggpubr` (figures); `vegan`, `dendextend`, `factoextra`, `mclust`
(multivariate / clustering analyses).

## License

MIT.
