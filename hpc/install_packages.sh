#!/bin/bash
# Installs the R packages of the pipeline into R_LIBS_USER, once, on the login node (the
# compute nodes may have no internet), from the repo root:
#   bash hpc/install_packages.sh [ppindices source folder, default /work/haarscheid/ppindices]
# Keep R_MODULE and R_LIBS_USER the same as in hpc/resolution_gwas.slurm.
set -euo pipefail
R_MODULE="lang/R/4.3.3-gfbf-2023b"
export R_LIBS_USER="/work/haarscheid/R/plastquant-4.3"
PPINDICES_SRC="${1:-/work/haarscheid/ppindices}"
module load "$R_MODULE"
mkdir -p "$R_LIBS_USER"
Rscript --vanilla - <<RSCRIPT
lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib, .libPaths()))
cran <- c("here", "rrBLUP", "ggplot2", "pheatmap", "rpart", "agricolae", "emmeans", "VIM")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
install.packages("$PPINDICES_SRC", lib = lib, repos = NULL, type = "source")
ok <- vapply(c(cran, "ppindices"), requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("not installed: ", paste(names(ok)[!ok], collapse = ", "))
stopifnot(ppindices::calculate_reaction_norm_slope(c(8, 6, 4, 2)) > 0)
cat("All packages installed in", lib, "\n")
RSCRIPT
