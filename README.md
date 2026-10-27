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
  variance/                latent variance components (target and realised h2, share of the causal SNPs)
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
and `variance/`.

All settings of a simulation are arguments of `sim_config()`,
defined in one file, `R/sim/simulation_config.R. 
The same file holds the shapes and their parameters (`shape_parameters()`) and the
ranges and distributions of the parameters (`default_study_space()`); to change a range,
edit `default_study_space()` or pass a modified copy as `study_space`. A few fixed model
choices are in the model code itself: allele frequencies U(0.2, 0.8) and effects ±1 in
`R/sim/snp_genotypes.R`, the shape functions in `R/sim/reaction_norm_shapes.R`.

Every argument with its default, to copy and adapt:

```r
cfg <- sim_config(
  h2           = 0.5,                    
  Mc           = 5,                      
  q            = 0.5,
  stability_sd = 0.2,
  n_per_shape  = 200,
  n_rep        = 3,
  n_snps       = 200000,
  ld_strength  = 0,
  ld_r2        = 0.1,
  n_points     = 6,
  seed         = 1,
  sigma_L      = 0.1,
  error_share  = 0.02,
  study_space  = default_study_space()
)
```

The settings:

Required:

- `h2`: narrow-sense heritability of every reaction norm parameter on the latent scale
  (see *Heritability* below)
- `Mc`: causal SNPs per parameter

With defaults:

- `q` (default 0.5): share of the additive genetic variance that comes from the causal
  SNPs; the rest is polygenic background. Each causal SNP explains `q · h2 / Mc` of the
  variance
- `stability_sd` (default 0.2): how much the genotypes differ in the spread of their
  replicates, as the SD of ln ν_i over the genotypes (ν_i: within-genotype variance of
  genotype i relative to the average, mean 1). 
- `n_per_shape` (default 200): genotypes per shape
- `n_rep` (default 3): replicate plants per genotype
- `n_snps` (default 200000): simulated SNPs
- `ld_strength` (default 0): linkage disequilibrium (LD) between neighbouring SNPs. With 0 every SNP is drawn
  independently (no LD). With LD, the causal SNPs of a shape lie more than two LD windows
  apart. `R/analysis/ld_decay.R` plots the mean r² against the distance
  (`plots/ld_decay.pdf`)
- `ld_r2` (default 0.1): sets the LD window, the largest distance in SNPs at which the
  mean r² is still at least `ld_r2`; measured from the simulated SNPs and stored in
  `data/config.rds` (`ld_window`). Only used with `ld_strength` > 0
- `n_points` (default 6): measurement points used by `simulate_reaction_norms()` when
  working in memory; `sample_resolution()` takes its own number of points
- `seed` (default 1): the same seed gives the same genetics for every `h2`
- `sigma_L` (default 0.1): SD of ln L over the plants
- `error_share` (default 0.02): measurement error as share of the trait variance
- `study_space` (default `default_study_space()`): range and distribution of every
  bounded parameter

**Heritability.** Each reaction norm parameter θ (L, a, c, w, k, ψ) is generated on a
latent scale with phenotypic variance 1 (`R/sim/genetic_values.R`). For genotype i and
replicate plant r:

| Quantity | Definition | Variance |
|---|---|---|
| genetic value | G_i = √h2 · (√q · S_i + √(1−q) · B_i) | V_A = h2 |
| causal SNPs | √(h2·q) · S_i (S: SNP score of the causal SNPs) | V_Q = q · h2 |
| polygenic background | √(h2·(1−q)) · B_i | (1 − q) · h2 |
| stability of genotype i | σ_i² = exp(s · ℓ_iV − s²/2), s = `stability_sd` | mean 1 |
| plant | ℓ_ir = G_i + √((1−h2) · σ_i²) · ζ_ir | V_P = 1 |

S, B and ζ are standard normal and independent; S and B are additive, so V_A = V_G. The
within-genotype variance is V_W = (1 − h2) · E[σ_i²] = 1 − h2, so **h2 = V_A / V_P is the
narrow-sense heritability on the latent scale**, exact in the population (in a finite
population the realised values deviate by chance; `variance/variance_latent.csv` lists
target and realised values). The stability ℓ_iV has its own causal SNPs (V in the truth
table); they change only the spread of the replicates, not their mean, so they add
nothing to V_A. The repeatability of genotype i is ρ_i = h2 / (h2 + (1 − h2) · σ_i²)
(column `rho` of the parameter tables). The parameters are then mapped from the latent
scale to their ranges (`study_space`); on the parameter scale and for the PPIs the
heritability is therefore no longer exactly h2.

### 1b. Alternatively: simulate on a real SNP panel

Simulated SNPs have idealised LD. For realistic LD, allele frequencies and relatedness, the
reaction norms can instead be simulated on the SNPs of real lines (for example the maize
panel in `data/A-Datasets/`):

```r
simulate_population_from_snps(sim_config(h2 = 0.4, Mc = 10), snps, POPULATION_DIR)
```

`snps` is a matrix (lines × SNPs, dosages 0/1/2, line IDs as row names, SNP IDs as column
names, no missing values) or a file with it (`.rds`, or `.csv` with a first column
`Genotype`). The result is a population folder like the one of step 1, and all following
steps run unchanged. What differs from step 1 (`R/sim/snp_panel.R`):

- **Genotypes:** every shape is simulated on all lines of the panel. The shapes stay
  separate (their own reaction norms, PPIs and GWAS); `n_per_shape` and `n_snps` of
  `sim_config()` are replaced by the panel.
- **Causal SNPs:** only SNPs with minor allele frequency ≥ `min_maf` (0.05), and within a
  shape no two causal SNPs in LD (r² ≥ `max_r2`, 0.1). L, a and V have the same causal SNPs
  and effects in every shape, so the shapes can be compared; the shape parameters (c, w,
  k, psi) get their own.
- **Standardisation:** the dosages of a causal SNP are centred and scaled with their own
  mean and SD over the lines. Simulated SNPs use 2p and √(2p(1−p)), which assumes
  Hardy–Weinberg; real lines are often inbred (almost only 0 and 2), so their variance
  would be wrong.
- **Hits in LD:** for every causal SNP, `data/causal_tags.csv` lists the SNPs with
  r² ≥ `max_r2` to it (its tag SNPs). In the GWAS evaluation a significant tag SNP finds
  the causal SNP. False positives are the hits that tag no causal SNP of the shape, counted
  as loci: hits with r² ≥ `max_r2` to each other count once.


**Example: the maize panel shipped in `data/A-Datasets/`** (403 lines; here only the
first 5000 SNPs, for a quick run):

```r
geno <- new.env()
load(here::here("data", "A-Datasets", "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)
snps <- geno$mat2[, 1:5000]   # lines x SNPs, dosages 0/1/2; geno$mat2 for all 978134 SNPs
rm(geno)

dir <- here::here("output", "populations", "maize_panel_test")
simulate_population_from_snps(sim_config(h2 = 0.4, Mc = 10), snps, dir)   # -> data/, variance/
res <- sample_resolution(dir, 6)                                          # -> n_6/data/
compute_population_ppis(res, "mean_before")                               # -> n_6/mean_before/ppis/
```

The GWAS follows as in step 4. With many SNPs it takes long: run it from the terminal on
several cores (see *Large populations* in step 4).

**Fewer SNPs: tag SNPs.** The full maize panel (978 134 SNPs) needs more memory than a
laptop has. `prune_snps()` reduces it to a subset of practically independent SNPs (tag
SNPs): along each chromosome, a SNP is kept only if its r² with every SNP already kept
among the preceding `window` SNPs (default 50) is below `r2` (default 0.1); otherwise a
kept SNP stands for it. `chrom` gives the chromosome of every SNP (default: all on one). The
causal SNPs are then drawn among the tag SNPs. Allele frequencies and relatedness of the
lines stay real; the LD between the markers is largely removed, so a causal SNP is found
by itself rather than together with its neighbours. The SNPs must be in genomic order:

```r
geno <- new.env()
load(here::here("data", "A-Datasets", "Genotyping_matrix_Negro_et_al_2019.Rdata"), envir = geno)
map <- readRDS(here::here("data", "A-Datasets", "SNPs_DL_windows_Negro_et_al_2019.rds"))
map <- map[map$SNP.name %in% colnames(geno$mat2), ]
map <- map[order(map$Chromosome, map$Position), ]            # genomic order
snps <- prune_snps(geno$mat2[, map$SNP.name], r2 = 0.1, window = 50, chrom = map$Chromosome)
rm(geno, map)
ncol(snps)                                                  # number of tag SNPs
```

`snps` then goes into `simulate_population_from_snps()` as above.

### 2. Measure it at a resolution

```r
RESOLUTION_DIR <- sample_resolution(POPULATION_DIR, 6)    # -> n_6/
```

The true norms are measured at 6 equidistant environments (1 to 10), and the
measurement error is added (2 % of the trait variance at these points). The same plants
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
With LD (`ld_strength` > 0), a causal SNP counts as found if a significant SNP lies within
its LD window, and the false positives are counted as loci: significant SNPs outside the
windows of all causal SNPs of the shape, where SNPs at most one window apart form one locus.

For the other handling, set `REPLICATE_HANDLING <- "replicates"` and run steps 3 and 4
again; the results go to `n_6/replicates/`.

**Large populations.** With many genotypes or SNPs the GWAS needs a lot of memory and
time, and an RStudio session can crash. Run it from the terminal with `Rscript` instead,
and on several cores with `n_cores` (the SNP tests are split over forked processes, which
RStudio does not handle well; the results are the same as on one core):

```bash
cd ~/PlastQuant
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla -e '
source(here::here("R", "pipeline", "load_pipeline.R"))
dir <- "output/populations/My_population/n_6"
compute_population_ppis(dir, "mean_before")            # -> n_6/mean_before/ppis/
run_population_gwas(dir, "mean_before", n_cores = 8)    # -> n_6/mean_before/gwas/, plots/
'
```

- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`: one thread per process, so the processes do
  not compete for the cores.
- `--vanilla`: R starts without `.Rprofile` and saved workspace; the packages are taken
  from the normal R library.
- `n_cores` (default 1): at most the number of cores of the machine, minus one or two for
  the system.

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
Mc ∈ {5, 10, 20} (narrow-sense heritability, `q` = 0.5) and run steps 1–4 at 6 points each:

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
