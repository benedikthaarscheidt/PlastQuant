# =============================================================================
# plotting_resolution.R — resolution correlation plots
# =============================================================================
# WHAT IT DOES: Reads the per-resolution correlation CSVs and renders the resolution
#   correlation plots across the defined reaction-norm ranges.
# REQUIRES:     correlations_*.csv files under the input directory (produced upstream).
# PRODUCES:     Resolution correlation figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plotting_resolution.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   dir_path     path   folder holding the correlations_*.csv inputs                       [COMMON]
#   range_names  character vector, default c("full","partial1","partial2")  ranges plotted
# =============================================================================

library(dplyr)
library(ggplot2)

ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}
dir_path <- path.expand(here::here("output", "correlation_summary_stats"))
output_dir <- path.expand(here::here("output", "plots"))
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
          print(p)
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
        print(p)
      }
    }
    
    dev.off()
    message("Saved PDF: ", pdf_filename)
  }
}
