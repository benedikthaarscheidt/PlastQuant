# Figures

What the main figures of the PPI comparison show, how they are computed and which script
produces them. All figures are produced by `sampling_analysis_publication.sh` (repository
root) for the population `output/populations/ppi_comparison_maize200`: 200 maize lines per
reaction norm shape (linear, Gaussian, sinusoidal, wave), 3 replicates, simulated with
h² = 0.4, M_c = 5, q = 0.5 (see the README for the simulation). Unless stated otherwise,
the PPIs are those of the replicate-averaged reaction norms measured at 6 points
(`n_6/mean_before/`).

The figures have no captions. `options(plastquant.captions = TRUE)` before sourcing a
script (or `CAPTIONS=TRUE` in `sampling_analysis_publication.sh`) adds a short caption
under every figure.

| Figure | Script | File |
|---|---|---|
| PPI distributions, clustering, Mantel tests | `R/analysis/plasticity_score_realisation_plot.R` | `n_6/mean_before/plots/figure2_2.pdf` |
| Top-k overlap | `R/analysis/rank_corr_biplot.R` | `n_6/mean_before/plots/plasticity_comprehensive_analysis/comprehensive_plasticity_analysis.pdf` (page 17, stacked over k) |
| Sampling resolution | `R/analysis/resolution_analysis.R`, `R/analysis/combined_resolution_figure.R` | `plots/combined_resolution_figure.pdf` |
| Maize phenotypes | step 4 of `sampling_analysis_publication.sh` | `output/populations/maize_traits_figure/n_8/mean_before/plots/figure2_2_maize.pdf` |

Paths are relative to the population folder, except the maize figure.

## Distributions, clustering and cross-shape comparison of the PPIs

**(a)** Distribution of each of the 28 PPIs over all 800 reaction norms (boxplots). For
display, values above the 95th percentile of each PPI are omitted. The PPIs differ by
orders of magnitude in scale, so raw values cannot be compared between PPIs.

**(b)** Hierarchical clustering of the PPIs over all shapes pooled. Two PPIs are close if
they rank the genotypes alike: the distance is 1 − τ, with τ the Kendall rank correlation
of the two PPIs over the genotypes; complete linkage. The number of clusters k ∈ {3, …, 10}
is the one with the largest average silhouette width, averaged over the four shapes;
label colours mark the clusters. Symbols mark the functional group of each PPI, assigned
by its formula independently of the clustering: dispersion (×), slope (□), range (●).
PPIs joined low in the tree rank the genotypes almost identically; clusters that mix
symbols show that the formula does not predict which PPIs agree.

**(c)** Mantel correlation between the PPI distance matrices of every pair of shapes and of
each shape with all shapes pooled (combined); Kendall correlation, 999 permutations;
\* p < 0.05, \*\* p < 0.01, \*\*\* p < 0.001. High values: the PPIs relate to each other in
the same way whatever the shape of the reaction norms.

PPIs that cannot be computed for every genotype (e.g. RN and RNN at very few points) are
left out.

## PPIs and summary statistics of the reaction norms

`n_6/mean_before/plots/ppi_summary_statistics.pdf` (same script; table:
`n_6/mean_before/analyses/ppi_summary_statistics.csv`). Spearman correlation over all
genotypes (all shapes pooled) between each PPI and five summary statistics of the
replicate-averaged reaction norm: mean (level), min, max, range (size of the response) and
the slope of the least-squares regression on the environment values (direction). Colours
and Benjamini–Hochberg adjusted significance stars only; rows in the order of the pooled
dendrogram, labels in its cluster colours. It shows which property of the reaction norms
each cluster follows.

## PPIs grouped by their correlation profiles

`n_6/mean_before/plots/ppi_profile_clustering.pdf` (`R/analysis/parameter_ppi_correlation.R`,
functions in `R/analysis/ppi_profiles.R`; tables `analyses/ppi_profile_groups.csv` and
`analyses/ppi_profile_comparison.csv`). Every PPI is described by a profile, the vector of
its signed Spearman correlations (a) with the true parameters of every shape and (b) with the
five summary statistics above. The PPIs are clustered on each profile (Euclidean distance,
average linkage, number of groups by silhouette width); the name colours show these groups,
the squares the groups of the Kendall clustering of the PPIs (all shapes pooled). Three
comparisons, each by a Mantel test between the distance matrices (rank correlation of the
distances, 999 permutations) and the adjusted Rand index of the groups: (a) vs Kendall,
(b) vs Kendall and (a) vs (b). PPIs in one profile group measure the same properties of the
reaction norm; agreement with the Kendall groups shows that PPIs rank genotypes alike because
of that. Only for simulated populations (true parameters needed); for the maize data the
profile from the summary statistics alone is drawn by `plasticity_score_realisation_plot.R`
(`ppi_profile_clustering_summary_statistics.pdf`).

For the maize data (`maize_traits_figure/n_8/...`) the three traits are pooled. They differ
in units and scale, so correlations with the mean there partly reflect the differences
between the traits rather than between lines.

## Top-k overlap among the PPIs

Every PPI ranks the 800 reaction norms (all shapes pooled) by decreasing value; its top-k
set holds the k genotypes with the highest values. Each cell shows the share of genotypes
the top-k sets of two PPIs have in common, |A ∩ B| / k, from 0 (no shared genotypes) to 1
(identical sets). Rows and columns are ordered by average-linkage clustering of
1 − overlap. The PDF has one page per k ∈ {3, 5, 10, 15, 20, 30, 50, 75, 100, 150}
(pages 7–16; the title of each page gives its k).

How the overlap depends on k is shown on three further pages. All of them use the same
groups of PPIs: average-linkage clustering of 1 − the overlap averaged over all k, cut
into the number of groups (2–8) with the largest average silhouette width.

- **Page 17, stacked overlap matrices:** the matrices of k = 5, 15, 30, 75 and 150 as
  layers (smallest k at the bottom), all in the same PPI order, names along the front and
  the right edge. Blocks
  present in every layer are groups of PPIs that select the same genotypes whatever the
  number selected.
- **Page 18, distribution over the PPI pairs:** boxplots of the overlap of all PPI pairs at
  every k, split into pairs within a group and pairs between groups; dashed line: the
  overlap expected by chance, k / number of genotypes. Between-group boxes near the chance
  line mean the groups select unrelated genotypes; within-group boxes that stay high mean
  the groups hold at every k.
- **Page 19, per PPI:** the mean overlap of each PPI with all other PPIs at every k (in %).
  Dark rows are PPIs that select other genotypes than the rest; rows that brighten with k
  are PPIs that differ from the others only in the very top genotypes.

For the stability indices (PSI, RSI) a high value means low plasticity, so their top-k sets
hold the least plastic genotypes. Low overlap between two PPIs means that the choice
between them changes which genotypes are selected as most plastic.

## Sensitivity to the sampling resolution

The same plants are measured at n ∈ {2, 3, 4, 5, 6, 11, 26, 50} equidistant points
(folders `n_2` … `n_50`), so differences between resolutions come from the sampling alone.
At every resolution, the PPIs and nine summary statistics of each replicate-averaged
reaction norm are computed: minimum, maximum, mean, median, range, variance, the slope of
the least-squares regression on the environment values, and the mean of the lower and of
the upper half of the points. All shapes are pooled.

**(a)** For every pair of PPI and summary statistic, the Pearson correlation across
genotypes is computed at every resolution and regressed linearly on n. Tile colour: slope
of this regression (blue: the correlation decreases with denser sampling; red: it
increases); significance of the slope by its t-test, \* p < 0.05, \*\* p < 0.01,
\*\*\* p < 0.001. A non-zero slope means that the relation between the PPI and that property
of the reaction norm depends on the number of points.

**(b)** For every PPI, the linear mixed model PPI ~ n + nine summary statistics +
(1 | genotype), fitted over all resolutions (lmerTest), with the PPI z-standardised over
all genotypes and resolutions. Points: the conditional slope β̂ₙ with 95 % confidence
interval, i.e. the change of the PPI in standard deviations per additional sampling point
with the summary statistics and the genotype held fixed; being standardised, it is
comparable between PPIs. A PPI is sample-size-sensitive (red) if its Benjamini–Hochberg
adjusted p < 0.05 and |β̂ₙ| ≥ 0.001, otherwise robust (blue).

The tables behind the figure: `analyses/resolution/correlations.csv`,
`correlation_regression.csv` and `mixed_model.csv` in the population folder.

## PPIs of the maize phenotypes

The same figure as the first one, for real data: biomass, leaf area and water-use
efficiency (WUE) of 244 maize lines in eight environments (seasons S12, S13, S16, W13,
each under drought and well-watered conditions). The environments are not on a controlled
gradient; each gets an environment index as in Finlay–Wilkinson regression, the mean of
the trait over all lines in that environment, separately for every trait. The three traits
take the place of the shapes: (a) distributions over all lines and traits, (b) clustering
over all traits pooled, with k chosen over the three traits, (c) Mantel correlations
between the traits.

The data are loaded by `R/gwas/scenarios/scenario_maize.R`; step 4 of
`sampling_analysis_publication.sh` runs that script up to its PPI step only (no GWAS) and
writes into its own folder, so `output/populations/maize` is not touched.
