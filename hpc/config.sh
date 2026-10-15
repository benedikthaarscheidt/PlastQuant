# Settings of the resolution GWAS on the cluster, read by the job scripts in hpc/.
# Edit here, then submit with: bash hpc/submit_resolution_gwas.sh   (from the repo root)

# R from the module system (module avail lang/R) and the library for the packages
R_MODULE="lang/R/4.3.3-gfbf-2023b"
export R_LIBS_USER="$HOME/R/plastquant-4.3"

# Population: heritability, causal SNPs per parameter, seed (one population per run)
H2=0.6
MC=5
SEED=1

# Resolutions (measurement points) and replicate handling of the PPIs
RESOLUTIONS="2 3 4 5 6 11 26 50"
HANDLING="mean_before"            # or "replicates"

# Folder of the ppindices source (for hpc/install_packages.sh)
PPINDICES_SRC="$HOME/ppindices"
