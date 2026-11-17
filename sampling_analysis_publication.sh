#!/bin/bash
# sampling_analysis_publication.sh — PPI comparison and sampling-resolution analysis
# on 200 randomly drawn maize lines (tag SNPs), for the publication.
# Run from anywhere:  bash sampling_analysis_publication.sh
#   1. draw the lines and simulate the population on their maize genotypes
#   2. resolutions n = 2 ... 50 with PPIs and summary statistics (in parallel)
#   3. resolution analyses and the analyses at n = 6 (in parallel)
#   4. the PPIs of the real maize traits and their figure (supplementary figure S1)
# Logs per step in logs/sampling_analysis/.
set -uo pipefail
cd ~/PlastQuant

# ---- Settings ---------------------------------------------------------------
POP="output/populations/ppi_comparison_maize200"
SNPS="output/maize_tag_snps.rds"   # 403 lines x tag SNPs (prune_snps())
N_LINES=200                        # number of lines drawn at random from the panel
LINE_SEED=2026                     # seed of the random draw of the lines
NCORE=6                            # parallel R processes
MAIZE_POP="output/populations/maize_traits_figure"   # own folder for the real maize PPIs
CAPTIONS=FALSE                     # TRUE: short captions under the figures
# ------------------------------------------------------------------------------

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=logs/sampling_analysis
mkdir -p "$LOG"
say() { echo "$(date '+%F %T')  $*"; }

# 1. population on N_LINES randomly drawn maize lines
say "simulating $POP"
Rscript --vanilla -e "
source(here::here('R', 'pipeline', 'load_pipeline.R'))
snps <- readRDS('$SNPS')
set.seed($LINE_SEED)
lines <- sort(sample(rownames(snps), $N_LINES))
snps <- snps[lines, ]
simulate_population_from_snps(sim_config(h2 = 0.4, Mc = 5, stability_sd = 0.2, n_rep = 3,
                                         error_share = 0.02, seed = 1),
                              snps, '$POP')
writeLines(lines, file.path('$POP', 'data', 'maize_lines.txt'))
" > "$LOG/01_simulate.log" 2>&1 || { say "simulation FAILED (see $LOG/01_simulate.log)"; exit 1; }

# 2. resolutions in parallel: every process writes only its own n_<n> folder
say "resolutions on $NCORE cores"
printf '%s\n' 2 3 4 5 6 11 26 50 | xargs -P "$NCORE" -I{} sh -c '
  Rscript --vanilla -e "source(here::here(\"R\", \"pipeline\", \"load_pipeline.R\")); invisible(compute_resolutions(\"$1\", $2))" \
    > "$3/02_n$2.log" 2>&1 && echo "  n = $2 done" || echo "  n = $2 FAILED (see $3/02_n$2.log)"' _ "$POP" {} "$LOG"

# 3. analyses: independent groups in parallel, each in its own R process
run() {
  local name=$1; shift
  local code="options(plastquant.captions = $CAPTIONS)
POPULATION_DIR <- '$POP'; RESOLUTION_DIR <- file.path(POPULATION_DIR, 'n_6'); REPLICATE_HANDLING <- 'mean_before'"
  for s in "$@"; do code="$code
source(here::here('R', 'analysis', '$s.R'))"; done
  if Rscript --vanilla -e "$code" > "$LOG/03_$name.log" 2>&1
  then say "  $name done"; else say "  $name FAILED (see $LOG/03_$name.log)"; fi
}

say "analyses"
run resolution        resolution_analysis combined_resolution_figure correlation_resolution correlation_resolution_experiment &
run realisation       plasticity_score_realisation_plot &
run rank_agreement    rank_corr_biplot &
run sumstat_corr      correlation_summary_stats_PS &
run kendall           kendall_rank_corr kendall_tau_corr_clustering &
run reaction_norms    plotting_reaction_norms_all &
run parameters        parameter_ppi_correlation &
run clustering        clustering_rn_summary_stats &
run regression        regression_summarystats &
wait

# 4. real maize phenotypes (biomass, leaf area, WUE; 8 environments): PPIs and their figure.
#    R/gwas/scenarios/scenario_maize.R is read as it is and run only up to its PPI/GWAS
#    step (loading phenotypes, environments and SNPs, writing the population folder); the
#    PPIs are computed here, without GWAS, into an own folder (MAIZE_POP), so neither the
#    scenario script nor output/populations/maize is touched. Only the first 1000 SNPs are
#    stored, as the figure does not use them.
say "maize phenotypes"
Rscript --vanilla -e "
options(plastquant.captions = $CAPTIONS)
POPULATION_DIR <- here::here('$MAIZE_POP'); SNP_SUBSET <- 1000
code <- readLines(here::here('R', 'gwas', 'scenarios', 'scenario_maize.R'))
cut <- grep('compute_population_ppis(RESOLUTION_DIR)', code, fixed = TRUE)[1]
if (is.na(cut)) stop('scenario_maize.R: PPI step not found')
eval(parse(text = code[seq_len(cut - 1)]))              # phenotypes, SNPs, population folder
compute_population_ppis(RESOLUTION_DIR, 'mean_before')  # -> n_8/mean_before/ppis/
form_ranges <- list(biomass = NULL, leafArea = NULL, WUE = NULL)   # maize mode of the figure script
REPLICATE_HANDLING <- 'mean_before'
source(here::here('R', 'analysis', 'plasticity_score_realisation_plot.R'))
" > "$LOG/04_maize.log" 2>&1 && say "  maize done" || say "  maize FAILED (see $LOG/04_maize.log)"

say "done: results in $POP and $MAIZE_POP"
