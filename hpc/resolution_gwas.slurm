#!/bin/bash
#SBATCH --job-name=resolution_gwas
#SBATCH --partition=all,gpu      # on the partition "all" or "gpu" (whichever starts first)
#SBATCH --nodes=1                # on a single node
#SBATCH --ntasks=1               # with a single task
#SBATCH --cpus-per-task=8        # one CPU per resolution: the GWAS of the resolutions run in parallel
#SBATCH --mem=128G               # all parallel GWAS together (each holds the SNP matrix)
#SBATCH --time=1-00:00:00        # maximum job runtime in "d-hh:mm:ss"
#SBATCH --output=resolution_gwas-%j.out   # the standard output (and stderr) of the job

# Influence of the resolution on the GWAS, for one population. Submit from the repo root:
#   sbatch hpc/resolution_gwas.slurm
# Steps:
#   1. simulate the population (kept if it exists)
#   2. regression_sumstats_dataprep.R: all resolutions n_<n> with their PPIs and summary statistics
#   3. GWAS of every resolution, in parallel (one R process per CPU); log per resolution in logs/
#   4. combine_resolution_gwas.R: power over the resolutions -> <population>/analyses/resolution_gwas/

# ---- Settings ----------------------------------------------------------------
H2=0.6                  # share of the latent genetic variance from the causal SNPs
MC=5                    # causal SNPs per parameter
SEED=1                  # seed of the population
HANDLING="mean_before"  # replicate handling of the PPIs: "mean_before" or "replicates"
R_MODULE="lang/R/4.3.3-gfbf-2023b"
export R_LIBS_USER="$HOME/R/plastquant-4.3"   # packages from hpc/install_packages.sh
# ------------------------------------------------------------------------------

set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load "$R_MODULE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per R process: the CPUs go to the resolutions
NCPU=${SLURM_CPUS_PER_TASK:-1}
POP="output/populations/resolution_$(awk "BEGIN{printf \"%.0f\", $H2*100}")_${MC}_seed${SEED}"
mkdir -p logs
echo "$(date '+%F %T')  population $POP, $NCPU CPUs"

# 1. population
Rscript --vanilla R/gwas/scenarios/resolution_gwas.R simulate "$POP" "$H2" "$MC" "$SEED" || exit 1

# 2. resolutions with PPIs and summary statistics
Rscript --vanilla -e "POPULATION_DIR <- '$POP'; source(here::here('R', 'analysis', 'regression_sumstats_dataprep.R'))" || exit 1
echo "$(date '+%F %T')  resolutions ready"

# 3. GWAS of every resolution, at most NCPU at a time
ls -d "$POP"/n_* | xargs -P "$NCPU" -I{} sh -c '
  N=$(basename "$1")
  if Rscript --vanilla R/gwas/scenarios/resolution_gwas.R gwas "$1" "$2" > "logs/gwas_${N}-$3.log" 2>&1
  then echo "$(date "+%F %T")  GWAS $N done"
  else echo "$(date "+%F %T")  GWAS $N FAILED (see logs/gwas_${N}-$3.log)"; fi' _ {} "$HANDLING" "${SLURM_JOB_ID:-local}"

# 4. power over the resolutions (the failed ones are left out)
Rscript --vanilla R/gwas/scenarios/combine_resolution_gwas.R "$POP" "$HANDLING"
echo "$(date '+%F %T')  done"
