#!/bin/bash
# Installs the R packages of the pipeline into R_LIBS_USER (hpc/config.sh), once, on the
# login node, from the repo root:  bash hpc/install_packages.sh
# ppindices is installed from its source folder PPINDICES_SRC (copy it to the cluster first).
set -euo pipefail
source hpc/config.sh
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
