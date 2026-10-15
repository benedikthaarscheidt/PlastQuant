#!/bin/bash
#SBATCH --job-name=rg_simulate
#SBATCH --partition=all          # on the partition "all"
#SBATCH --nodes=1                # on a single node
#SBATCH --ntasks=1               # with a single task
#SBATCH --cpus-per-task=1        # with that many CPU-cores
#SBATCH --mem=40G                # and this amount of memory
#SBATCH --time=0-02:00:00        # maximum job runtime in "d-hh:mm:ss"
#SBATCH --output=logs/%x-%j.out  # standard output and error

# Step 1: simulate the population (once); run from the repo root
cd "$SLURM_SUBMIT_DIR"
source hpc/config.sh
module load "$R_MODULE"
Rscript --vanilla R/gwas/scenarios/resolution_gwas.R simulate "$H2" "$MC" "$SEED"
