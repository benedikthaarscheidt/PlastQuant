#!/bin/bash
#SBATCH --job-name=rg_combine
#SBATCH --partition=all          # on the partition "all"
#SBATCH --nodes=1                # on a single node
#SBATCH --ntasks=1               # with a single task
#SBATCH --cpus-per-task=1        # with that many CPU-cores
#SBATCH --mem=8G                 # and this amount of memory
#SBATCH --time=0-00:30:00        # maximum job runtime in "d-hh:mm:ss"
#SBATCH --output=logs/%x-%j.out  # standard output and error

# Step 3: power over all resolutions of the population; run from the repo root
cd "$SLURM_SUBMIT_DIR"
source hpc/config.sh
module load "$R_MODULE"
POP="output/populations/resolution_$(awk "BEGIN{printf \"%.0f\", $H2*100}")_${MC}_seed${SEED}"
Rscript --vanilla R/gwas/scenarios/combine_resolution_gwas.R "$POP" "$HANDLING"
