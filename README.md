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
                                  compute_population_ppis.R
  gwas/                           GWAS functions, run_population_gwas.R, power per shape,
                                  plot_gwas_power.R (power heatmaps), 06_analyze_gwas_results.R
                                  (scenario comparison), 06_analyze_maize_gwas_results.R
    scenarios/                    heritability scenarios, combine_scenarios.R, scenario_maize.R
                                  (the maize data through the pipeline)
  analysis/                       comparison of the PPIs: analyses and their figures,
                                  resolution stability, summary statistics, rank correlation
data/                             committed input datasets (maize phenotypes and genotypes)
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
computed from them. Below a resolution, the PPIs and everything computed from them sit
in one folder per replicate handling (see step 3), so the two never overwrite each other.

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
    analyses/  plots/        of the measurement itself (summary statistics, measured reaction norms)
    mean_before/             PPIs from the replicate-averaged reaction norms (default)
      ppis/  gwas/  analyses/  plots/
    replicates/              PPIs from every replicate's own reaction norm
      ppis/  gwas/  analyses/  plots/
  n_11/ ...                same population, other resolutions
  analyses/  plots/        analyses and figures across resolutions
```

## Running the pipeline

Load the functions once per R session, from the project root (again after the code
has changed):

```r
setwd("~/PlastQuant")
here::i_am("R/pipeline/load_pipeline.R")   # resets the root to ~/PlastQuant
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

### 3. Choose the replicate handling and compute the PPIs

The PPIs of a resolution can be computed in two ways. Choose one by setting
`REPLICATE_HANDLING`, and pass it to every following step; all results of that choice
go into its own folder below the resolution (`n_6/mean_before/` or `n_6/replicates/`):

- `"mean_before"` (the default): the replicates of a genotype are averaged per
  environment first, and the PPI is computed on this averaged reaction norm — one value
  per genotype.
- `"replicates"`: the PPI is computed on every replicate's own measured reaction norm,
  nothing is averaged — one value per replicate (in the GWAS, the replicates of a
  genotype are then repeated observations of the same line).

```r
REPLICATE_HANDLING <- "mean_before"                             # or "replicates"
compute_population_ppis(RESOLUTION_DIR, REPLICATE_HANDLING)     # all 28 PPIs -> n_6/mean_before/ppis/
```

Both handlings can be computed for the same resolution; neither overwrites the other.
Without the second argument, `"mean_before"` is used. FW uses the reaction norms of the
whole population and is always computed once on the averaged norms, one value per
genotype, with either handling.

The console shows, for every PPI and shape, the number of genotypes and the time it
took, the output of the PPI functions themselves (e.g. the progress bars of PPF) and
every error (with the PPI, the genotype and the message; that value becomes `NA`).

To leave PPIs out, use `exclude`; to add a PPI, register it and pass the registry:

```r
compute_population_ppis(RESOLUTION_DIR, REPLICATE_HANDLING, exclude = c("PPF"))
reg <- register_ppi(default_ppi_registry(),
                    ppi_entry("MyIndex", function(trait_values, env_values) diff(range(trait_values)) / diff(range(env_values)),
                              inputs = c(trait_values = "trait", env_values = "env")))
compute_population_ppis(RESOLUTION_DIR, REPLICATE_HANDLING, registry = reg)
```

### 4. GWAS (optional)

The GWAS uses the PPIs of the chosen handling (step 3 must have been run for it):

```r
run_population_gwas(RESOLUTION_DIR, REPLICATE_HANDLING)   # -> n_6/mean_before/gwas/, n_6/mean_before/plots/
```

rrBLUP mixed model with the kinship matrix, one GWAS per PPI and shape
(`gwas/gwas_results.csv`). If the population has a truth table (always for simulated
populations), the same step evaluates it: the power per shape (`gwas/gwas_power.csv`;
Bonferroni `alpha` / number of SNPs, `alpha = 0.05` unless given) and its heatmaps
(`plots/gwas_power_heatmap.pdf`). A parameter is evaluated within its own shape, and a hit
on a SNP that is not causal for that shape counts as a false positive.

For the other handling, set `REPLICATE_HANDLING <- "replicates"` and run steps 3 and 4
again; the results go to `n_6/replicates/`.

### Figures and analyses of one resolution

Set `RESOLUTION_DIR` (step 2) and source the script. The scripts that use PPIs read and
write `RESOLUTION_DIR/<REPLICATE_HANDLING>/`; `REPLICATE_HANDLING` is `"mean_before"` unless
you set it (the analyses need one value per genotype, so `"mean_before"`). The measured
reaction norms go to `RESOLUTION_DIR/plots/`. Steps 1–3 are needed.

```r
source(here::here("R", "analysis", "plasticity_score_realisation_plot.R"))  # Fig. 2: distributions, clustering, Mantel
source(here::here("R", "analysis", "rank_corr_biplot.R"))                   # Fig. 3: rank agreement, top-k overlap
source(here::here("R", "analysis", "plotting_reaction_norms_all.R"))        # measured reaction norms and histograms
source(here::here("R", "analysis", "correlation_summary_stats_PS.R"))       # PPIs vs summary statistics
source(here::here("R", "analysis", "regression_summarystats.R"))           # PPIs regressed on summary statistics
source(here::here("R", "analysis", "kendall_rank_corr.R"))                 # rank tables ...
source(here::here("R", "analysis", "kendall_tau_corr_clustering.R"))       # ... and their agreement analyses
source(here::here("R", "analysis", "clustering_rn_summary_stats.R"))       # clustering of the reaction norms
source(here::here("R", "analysis", "summary_stats_vs_scores.R"))           # functions relating PPIs to summary statistics
source(here::here("R", "analysis", "parameter_ppi_correlation.R")).        # correlation analysis between traits (PPI realisations) and ground truth reaction norm parameters
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

### The maize data

`scenario_maize.R` runs the empirical maize data (244 varieties; biomass, leaf area and
WUE in 8 environments) through the same steps:

```bash
Rscript --vanilla R/gwas/scenarios/scenario_maize.R       # -> output/populations/maize/
```

1. **Load and structure.** It reads the three traits and the genotypes of the varieties
   from `data/A-Datasets/` and writes them as a population folder: the SNPs (978,134,
   stored as `.rds`) in `maize/data/`, the measurements as the resolution `maize/n_8/`.
   The environments are not controlled, so each environment gets an environment index
   as in Finlay–Wilkinson: for each trait its own, the environment mean of that trait as
   given in the FW/AMMI files of the data (`ENV_Biomass`, `ENV_Leaf_area`, `ENV_WUE`).
2. **PPIs and GWAS.** `compute_population_ppis()` and `run_population_gwas()` on `n_8`
   (`n_8/mean_before/`). There is no truth, so there is no power.
3. **GWAS summary.** `06_analyze_maize_gwas_results.R`: hits per PPI and trait
   (Bonferroni and −log10 P > 3, > 5), which PPI finds which SNP, the allele frequency of
   the hit SNPs, the similarity of the PPIs.

The maize figures read the maize PPIs from `maize/n_8/mean_before/`:

```r
source(here::here("R", "analysis", "plasticity_score_realisation_plot_maize.R"))   # -> maize/n_8/mean_before/plots/
```

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

If the environment values differ per trait (e.g. an environment index of each trait),
`environments.csv` has one value column per trait, named like the trait files, instead
of `Value`:

```
Column, biomass, leafArea
dry,    194.4,   0.26
wet,    263.0,   0.32
```

Then run steps 3 and 4 and the figures of one resolution on it:

```r
RESOLUTION_DIR <- here::here("output", "populations", "my_data", "n_5")
compute_population_ppis(RESOLUTION_DIR)    # -> n_5/mean_before/ppis/
run_population_gwas(RESOLUTION_DIR)        # needs snp_matrix.csv; -> n_5/mean_before/gwas/
```

Notes:

- Own data is measured at its environments only, so it cannot be measured at other
  resolutions; the stability analysis across resolutions needs a simulated population.
- Without `truth.csv` the GWAS gives results but no power.
- Some PPIs need a minimum number of environments (RNN: more than 3). A PPI that cannot
  be computed gives `NA`; the error is printed with the PPI and the genotype, and the
  run continues.
- A large SNP matrix can be stored as `snp_matrix.rds` (a matrix with the genotype IDs as
  row names) instead of `snp_matrix.csv`.
- Format errors are reported with the file and the offending column or genotype.

| Script | Simulated | Own data |
|---|---|---|
| PPIs, GWAS (`compute_population_ppis`, `run_population_gwas`) | yes | yes (GWAS with SNPs) |
| Power and its heatmaps (part of `run_population_gwas`) | yes | with SNPs and truth |
| Fig. 2, Fig. 3, summary-statistic and Kendall analyses | yes | yes |
| Measured reaction norms (`plotting_reaction_norms_all.R`) | yes | yes (genotype means instead of mother norms) |
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
