# =============================================================================
# plotting_resolution.R — resolution correlation plots
# =============================================================================
# WHAT IT DOES: Reads the per-resolution correlation CSVs and renders the resolution
#   correlation plots across the defined reaction-norm ranges.
# REQUIRES:     correlations_*.csv of resolution.R in POPULATION_DIR/analyses/correlation_summary_stats/.
# PRODUCES:     plots_<range>_<analysis>_per_score.pdf in POPULATION_DIR/plots/.
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "ppi_comparison")
#               source(here::here("R", "analysis", "plotting_resolution.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   POPULATION_DIR  path, REQUIRED (set before sourcing)   population folder                [COMMON]
#   range_names  character vector, default c("full","partial1","partial2")  ranges plotted
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

library(dplyr)
library(ggplot2)

ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}
if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")
source(here::here("R", "pipeline", "load_pipeline.R"))   # caption helpers
dir_path <- file.path(POPULATION_DIR, "analyses", "correlation_summary_stats")
output_dir <- file.path(POPULATION_DIR, "plots")
ensure_dir_exists(dir_path)
ensure_dir_exists(output_dir)


# Our predetermined range tags as they appear in your filenames.
# For example, files contain "full", "partial1", or "partial2".
range_names <- c("full", "partial1", "partial2")

# List all CSV files in the directory that follow our naming pattern.
files <- list.files(path = dir_path, pattern = "^correlations_.*\\.csv$", full.names = TRUE)

# Read and combine CSV files.
# For each file, extract the range name and the indices part from the filename.
data_all <- bind_rows(lapply(files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  base <- basename(f)
  # Extract range name (matches anything between "correlations_" and "_interval_")
  rn <- sub("correlations_(.*?)_interval_.*\\.csv", "\\1", base)
  df$range_name <- rn
  # Extract the indices substring (assumes filename has _indices_<numbers separated by underscore>.csv)
  indices_str <- sub(".*_indices_(.*)\\.csv", "\\1", base)
  indices_vec <- strsplit(indices_str, "_")[[1]]
  n_selected <- length(indices_vec)
  df$n_selected <- n_selected
  return(df)
}))

# Convert relevant columns to numeric.
data_all <- data_all %>%
  mutate(n_selected   = as.numeric(n_selected),
         correlation  = as.numeric(correlation),
         p_value      = as.numeric(p_value))

# Filter to keep only rows with significant correlations (p_value < 0.05)
data_all <- data_all %>% filter(p_value < 0.05)

# Get unique analysis types and score names.
analysis_types <- unique(data_all$analysis)
score_names    <- unique(data_all$score)

# Captions
range_txt <- c(full = "the whole gradient (grid points 1-50 of 50)",
               partial1 = "the first half of the gradient (grid points 1-25 of 50)",
               partial2 = "the second half of the gradient (grid points 25-50 of 50)")
res_setup <- function(rn, analysis) paste0(
  "Population ", basename(normalizePath(POPULATION_DIR)), "; reaction norms sampled within ", range_txt[[rn]],
  " at the numbers of points on the x-axis; PPIs and summary statistics of the replicate-averaged reaction norm of each genotype, ",
  if (analysis == "global") "all shapes pooled. " else "one shape. ",
  "Kendall tau between the PPI and each summary statistic over the genotypes; only significant correlations shown (permutation test, 999 permutations, p < 0.05).")
res_reading <- "Flat lines: the PPI follows the same reaction norm properties at every sampling density; changes along x: what the PPI measures depends on how many points are sampled. Missing points: not significant at that density."

# Loop over each range.
for (rn in range_names) {
  # Subset data for the current range.
  data_range <- data_all %>% filter(range_name == rn)
  
  # For each analysis type (global or local), create a separate PDF.
  for (anal in analysis_types) {
    data_sub <- data_range %>% filter(analysis == anal)
    
    # Name the PDF file with the range tag and analysis type.
    pdf_filename <- file.path(output_dir, paste0("plots_", rn, "_", anal, "_per_score.pdf"))
    pdf(pdf_filename, width = 8, height = 6)
    
    # For each score, create one plot.
    for (sc in score_names) {
      df_sc <- data_sub %>% filter(score == sc)
      if (nrow(df_sc) == 0) next  # Skip if no data for this score.
      
      if (anal == "local") {
        # For local analysis, plot separate plots per genotype form.
        unique_gf <- unique(df_sc$genotype_form)
        for (gf in unique_gf) {
          df_sc_gf <- df_sc %>% filter(genotype_form == gf)
          p <- ggplot(df_sc_gf, aes(x = n_selected, y = correlation, color = stat, group = stat)) +
            geom_line(size = 1) +
            geom_point(size = 2) +
            scale_x_continuous(name = "Number of Selected Samples", 
                               breaks = seq(0, max(data_all$n_selected, na.rm = TRUE), by = 10)) +
            scale_y_continuous(limits = c(-1, 1)) +
            labs(title = paste("Score:", sc, "(Local,", rn, ", Genotype:", gf, ")"),
                 y = "Kendall Correlation",
                 color = "Summary Statistic") +
            theme_minimal()
          print(with_caption(p, res_setup(rn, "local"), res_reading))
        }
      } else {
        # For global analysis, plot one line per summary statistic.
        p <- ggplot(df_sc, aes(x = n_selected, y = correlation, color = stat, group = stat)) +
          geom_line(size = 1) +
          geom_point(size = 2) +
          scale_x_continuous(name = "Number of Selected Samples", 
                             breaks = seq(0, max(data_all$n_selected, na.rm = TRUE), by = 10)) +
          scale_y_continuous(limits = c(-1, 1)) +
          labs(title = paste("Score:", sc, "(Global,", rn, ")"),
               y = "Kendall Correlation",
               color = "Summary Statistic") +
          theme_minimal()
        print(with_caption(p, res_setup(rn, "global"), res_reading))
      }
    }
    
    dev.off()
    message("Saved PDF: ", pdf_filename)
  }
}
