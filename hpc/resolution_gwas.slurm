#!/bin/bash
#SBATCH --job-name=rg_resolution
#SBATCH --partition=all          # on the partition "all"
#SBATCH --nodes=1                # on a single node
#SBATCH --ntasks=1               # with a single task
#SBATCH --cpus-per-task=1        # with that many CPU-cores
#SBATCH --mem=40G                # and this amount of memory
#SBATCH --time=0-04:00:00        # maximum job runtime in "d-hh:mm:ss"
#SBATCH --output=logs/%x-%A_%a.out  # one file per array task

# Step 2: one array task per resolution: measure, PPIs, GWAS; run from the repo root.
# The array index picks the resolution from RESOLUTIONS in hpc/config.sh.
cd "$SLURM_SUBMIT_DIR"
source hpc/config.sh
read -r -a N <<< "$RESOLUTIONS"
n_points=${N[$SLURM_ARRAY_TASK_ID]}
module load "$R_MODULE"
Rscript --vanilla R/gwas/scenarios/resolution_gwas.R resolution "$H2" "$MC" "$SEED" "$n_points" "$HANDLING"
