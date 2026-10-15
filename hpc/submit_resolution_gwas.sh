#!/bin/bash
# Submits the resolution GWAS with the settings of hpc/config.sh, from the repo root:
#   bash hpc/submit_resolution_gwas.sh            simulate, all resolutions, combine
#   bash hpc/submit_resolution_gwas.sh --no-sim   reuse the simulated population
# The resolutions start when the simulation has finished successfully, the combination
# when all resolutions have finished (afterany: a failed resolution is left out).
set -euo pipefail
source hpc/config.sh
mkdir -p logs
read -r -a N <<< "$RESOLUTIONS"
last=$(( ${#N[@]} - 1 ))

if [[ "${1:-}" == "--no-sim" ]]; then
  res=$(sbatch --parsable --array=0-$last hpc/resolution_gwas.slurm)
else
  sim=$(sbatch --parsable hpc/simulate.slurm)
  echo "simulation:  job $sim"
  res=$(sbatch --parsable --dependency=afterok:$sim --array=0-$last hpc/resolution_gwas.slurm)
fi
echo "resolutions: job $res (n = $RESOLUTIONS)"
comb=$(sbatch --parsable --dependency=afterany:$res hpc/combine.slurm)
echo "combination: job $comb"
echo "follow with: squeue -u $USER;  logs in logs/"
