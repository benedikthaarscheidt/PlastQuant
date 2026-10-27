# =============================================================================
# summary_stats_vs_scores.R — relate summary statistics to plasticity scores
# =============================================================================
# WHAT IT DOES: Computes and visualises the relationship between reaction-norm summary
#   statistics and the plasticity scores (correlation matrices, distance measures,
#   clustering). Reads the PPIs and reaction norms of RESOLUTION_DIR.
# REQUIRES:     RESOLUTION_DIR with its PPIs (compute_population_ppis()).
# PRODUCES:     Summary-stat vs score figures under RESOLUTION_DIR/<REPLICATE_HANDLING>/analyses/summary_stats_vs_scores/.
# HOW TO RUN:   setwd("~/PlastQuant"); RESOLUTION_DIR <- here::here("output", "populations", "<name>", "n_6")
#               source(here::here("R","analysis","summary_stats_vs_scores.R"))
# -----------------------------------------------------------------------------
# PARAMETERS
#   RESOLUTION_DIR  path, REQUIRED (set before sourcing)   resolution folder               [COMMON]
#   REPLICATE_HANDLING  "mean_before" (default) or "replicates"   which PPIs of the resolution
#   Analysis settings are set in the body of the script.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

# PPIs and replicate-mean reaction norms of the resolution RESOLUTION_DIR
if (!exists("RESOLUTION_DIR")) stop("set RESOLUTION_DIR (a resolution folder, e.g. output/populations/<name>/n_6) first")
if (!exists("REPLICATE_HANDLING")) REPLICATE_HANDLING <- "mean_before"   # PPIs of mean_before/ or replicates/
PPI_SET_DIR <- file.path(RESOLUTION_DIR, REPLICATE_HANDLING)
source(here::here("R", "pipeline", "load_pipeline.R"))
analysis_inputs <- score_analysis_inputs(RESOLUTION_DIR, REPLICATE_HANDLING)

library(factoextra)
library(dbscan)
library(cluster)
library(mclust)
library(pheatmap)
library(corrplot)
library(vegan)
library(energy)
library(gridExtra)
library(grid)
library(gtable)

sig_stars <- function(p) {
  if (is.na(p)) {
    return("NA")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Define your scores_list
scores_list    <- analysis_inputs$scores_list      # one table per PPI: Score, Type (shape), Genotype_ID
rn_matrix_list <- analysis_inputs$rn_matrix_list   # per shape: genotypes x environments (replicate means)
genotype_forms <- analysis_inputs$genotype_forms



prepare_score_vector <- function(score_df, genotype_form) {
  # Assumes score_df has columns "Score" and "Type"
  subset(score_df, score_df[, "Type"] == genotype_form)[, "Score"]
}

# calc_similarity() and calc_similarity_perm() are in R/analysis/similarity.R (loaded by 03)

compute_summary_stats <- function(rn_matrix) {
  n <- ncol(rn_matrix)
  if(n < 2) stop("Reaction norm matrix must have at least 2 columns.")
  x_vals <- seq_len(n)
  slopes <- apply(rn_matrix, 1, function(x) coef(lm(x ~ x_vals))[2])
  ranges <- apply(rn_matrix, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  variances <- apply(rn_matrix, 1, function(x) var(x, na.rm = TRUE))
  
  half <- floor(n/2)
  mean_lower <- apply(rn_matrix, 1, function(x) mean(x[1:half], na.rm = TRUE))
  mean_upper <- apply(rn_matrix, 1, function(x) mean(x[(half+1):n], na.rm = TRUE))
  
  data.frame(
    genotype    = rownames(rn_matrix),
    min         = apply(rn_matrix, 1, min, na.rm = TRUE),
    max         = apply(rn_matrix, 1, max, na.rm = TRUE),
    mean        = apply(rn_matrix, 1, mean, na.rm = TRUE),
    median      = apply(rn_matrix, 1, median, na.rm = TRUE),
    slope       = slopes,
    range       = ranges,
    variance    = variances,
    mean_lower  = mean_lower,
    mean_upper  = mean_upper
  )
}


# Compute nine summary statistics for each row (genotype) in a reaction norm matrix.
compare_summary_stats_with_scores <- function(rn_matrix_list, scores_list, genotype_forms,
                                              global = FALSE,
                                              measures = c("kendall"), n_perm = 999,
                                              out_dir = file.path(PPI_SET_DIR, "analyses", "summary_stats_vs_scores")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Define the nine summary statistics.
  stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")
  score_names <- names(scores_list)

  # Captions
  setup <- paste(resolution_setup(RESOLUTION_DIR, REPLICATE_HANDLING,
                                  vapply(rn_matrix_list[genotype_forms], nrow, integer(1))),
                 "Summary statistics of the replicate-averaged measured reaction norm of each genotype (slope against the index of the measurement point).")
  perm_txt <- function(meth) paste0(meth, " correlation with a permutation test (", n_perm, " permutations; * p < 0.05, ** p < 0.01, *** p < 0.001).")
  scatter_reading <- "Shows which reaction norm properties a PPI follows and whether the relation is monotone, driven by a few genotypes or split by shape."
  table_reading <- "The marked statistics are the properties each PPI mainly measures; PPIs with the same marked rows measure the same properties."
  
  # Define a custom theme with reduced padding for a more filled cell.
  custom_theme <- ttheme_minimal(
    core = list(
      fg_params = list(col = "black", cex = 0.3),
      bg_params = list(fill = "grey95", col = "white"),
      padding = unit(c(2, 2), "mm")
    ),
    colhead = list(
      fg_params = list(col = "black", cex = 0.3, fontface = "bold"),
      bg_params = list(fill = "grey90", col = "white"),
      padding = unit(c(2, 2), "mm")
    ),
    rowhead = list(
      fg_params = list(col = "black", cex = 0.3, fontface = "italic"),
      bg_params = list(fill = "grey90", col = "white"),
      padding = unit(c(2, 2), "mm")
    )
  )
  
  if (!global) {
    message("Performing LOCAL (per genotype form) analysis...")
    pdf(file.path(out_dir, "scores_vs_summary_stats_local.pdf"), width = 10, height = 8)
    
    for (meth in measures) {
      for (form in genotype_forms) {
        message("Local Analysis - Method: ", meth, " | Form: ", form)
        rn_matrix <- rn_matrix_list[[form]]
        stat_df <- compute_summary_stats(rn_matrix)
        
        ## --- Scatter Plots ---
        for (score_name in score_names) {
          score_df <- scores_list[[score_name]]
          score_vector <- as.numeric(as.character(prepare_score_vector(score_df, form)))
          
          # Check if the score vector is empty or all NA
          if (length(score_vector) == 0 || all(is.na(score_vector))) {
            message("Skipping ", score_name, " for form ", form, ": score vector is empty or all NA.")
            next
          }
          
          # If the score vector has names, try to align it with stat_df$genotype.
          if (!is.null(names(score_vector))) {
            common_genos <- intersect(names(score_vector), stat_df$genotype)
            if (length(common_genos) == 0) {
              message("No matching genotypes for score ", score_name, " in form ", form, ". Skipping.")
              next
            }
            score_vector <- score_vector[common_genos]
            stat_df_sub <- stat_df[stat_df$genotype %in% common_genos, ]
          } else {
            stat_df_sub <- stat_df
          }
          
          # Check that the lengths match
          if (length(score_vector) != nrow(stat_df_sub)) {
            message("Length mismatch for score ", score_name, " in form ", form, 
                    ": score_vector length=", length(score_vector),
                    ", expected=", nrow(stat_df_sub), ". Skipping.")
            next
          }
          
          par(mfrow = c(3, 3)); caption_margin(4)
          for (st in stat_names) {
            x_vals <- stat_df_sub[[st]]
            y_vals <- score_vector
            test <- calc_similarity_perm(x_vals, y_vals, method = meth, n_perm = n_perm)
            corr_val <- test$estimate
            p_val <- test$p.value
            sig <- sig_stars(p_val)
            plot(x_vals, y_vals,
                 main = paste(form, ":", score_name, "vs.", st, "\n", meth, "corr =", round(corr_val, 3), "(", sig, ")"),
                 xlab = st, ylab = score_name,
                 pch = 19, col = "blue")
            legend("topleft", legend = paste(meth, "=", round(corr_val, 3), "(", sig, ")"), bty = "n")
          }
          base_caption(paste(setup, paste0("Shape: ", form, "."), "PPI", score_name, "against each summary statistic; each point a genotype;", perm_txt(meth)),
                       scatter_reading)
          par(mfrow = c(1, 1))
        }
        
        ## --- Correlation Matrix Table ---
        cor_mat <- matrix(NA_real_, nrow = length(stat_names), ncol = length(score_names),
                          dimnames = list(stat_names, score_names))
        pval_mat <- matrix(NA_real_, nrow = length(stat_names), ncol = length(score_names),
                           dimnames = list(stat_names, score_names))
        for (sc_idx in seq_along(score_names)) {
          sc_name <- score_names[sc_idx]
          score_df <- scores_list[[sc_name]]
          score_vector <- as.numeric(as.character(prepare_score_vector(score_df, form)))
          if (!is.null(names(score_vector))) {
            common_genos <- intersect(names(score_vector), stat_df$genotype)
            if (length(common_genos) == 0) {
              message("No matching genotypes for score ", sc_name, " in form ", form, ". Skipping matrix for this score.")
              next
            }
            score_vector <- score_vector[common_genos]
            stat_df_sub <- stat_df[stat_df$genotype %in% common_genos, ]
          } else {
            stat_df_sub <- stat_df
          }
          for (i in seq_along(stat_names)) {
            st <- stat_names[i]
            x_vals <- stat_df_sub[[st]]
            test <- calc_similarity_perm(score_vector, x_vals, method = meth, n_perm = n_perm)
            cor_mat[i, sc_name] <- round(test$estimate, 3)
            pval_mat[i, sc_name] <- round(test$p.value, 3)
          }
        }
        
        comb_mat <- matrix(NA_character_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                           dimnames = dimnames(cor_mat))
        for (i in seq_len(nrow(cor_mat))) {
          for (j in seq_len(ncol(cor_mat))) {
            comb_mat[i, j] <- paste0(cor_mat[i, j], " (", sig_stars(pval_mat[i, j]), ")")
          }
        }
        
        grid.newpage()
        grid.text(paste("Local Correlation Matrix - Method:", meth, "Form:", form),
                  gp = gpar(fontsize = 10, fontface = "bold"), y = 0.97)
        
        tbl <- tableGrob(comb_mat, rows = rownames(comb_mat), cols = colnames(comb_mat),
                         theme = custom_theme)
        
        
        for (j in seq_along(score_names)) {
          sc_name <- score_names[j]
          col_vals <- cor_mat[, sc_name]
          top_indices <- order(abs(col_vals), decreasing = TRUE)[1:2]
          highest_cell <- which(tbl$layout$name == "core-bg" & tbl$layout$t == (top_indices[1] + 1) & tbl$layout$l == (j + 1))
          second_cell  <- which(tbl$layout$name == "core-bg" & tbl$layout$t == (top_indices[2] + 1) & tbl$layout$l == (j + 1))
          
          if (length(highest_cell) > 0) {
            tbl$grobs[[highest_cell]]$gp$fill <- "green"
          }
          if (length(second_cell) > 0) {
            tbl$grobs[[second_cell]]$gp$fill <- "yellow"
          }
        }
        
        pushViewport(viewport(width = 0.9, height = 1))
        grid.draw(tbl)
        popViewport()
        grid_page_caption(paste(setup, paste0("Shape: ", form, "."), perm_txt(meth), "Rows: summary statistics, columns: PPIs; green: strongest, yellow: second strongest absolute correlation of each PPI."),
                          table_reading)
        
        for (sc_idx in seq_along(score_names)) {
          sc_name <- score_names[sc_idx]
          col_vals <- cor_mat[, sc_name]
          top_indices <- order(abs(col_vals), decreasing = TRUE)[1:2]
          message("LOCAL - Form: ", form, " | Method: ", meth, " | Score: ", sc_name,
                  " => Top correlations: ", paste(rownames(cor_mat)[top_indices], collapse = ", "))
        }
      }
    }
    dev.off()
    message("Local PDF saved as 'scores_vs_summary_stats_local.pdf'")
    
  } else {
    message("Performing GLOBAL (combined forms) analysis...")
    # Combine summary stats from all forms.
    combined_stats <- data.frame()
    for (form in genotype_forms) {
      rn_matrix <- rn_matrix_list[[form]]
      form_stats <- compute_summary_stats(rn_matrix)
      form_stats$Form <- form
      combined_stats <- rbind(combined_stats, form_stats)
    }
    # Define symbols for each genotype form.
    formSymbols <- c(linear = 16, gaussian = 17, sinusoidal = 18, wave = 15)
    
    pdf(file.path(out_dir, "scores_vs_summary_stats_global.pdf"), width = 12, height = 8)
    
    for (meth in measures) {
      message("Global Analysis - Method: ", meth)
      
      ## --- Scatter Plots ---
      for (score_name in score_names) {
        score_df <- scores_list[[score_name]]
        overall_scores <- c()
        for (form in genotype_forms) {
          score_subset <- as.numeric(as.character(prepare_score_vector(score_df, form)))
          overall_scores <- c(overall_scores, score_subset)
        }
        if(length(overall_scores) == 0 || all(is.na(overall_scores))) {
          message("Skipping ", score_name, " in global analysis: overall score vector is empty or all NA.")
          next
        }
        par(mfrow = c(3, 3)); caption_margin(4)
        for (st in stat_names) {
          x_vals <- combined_stats[[st]]
          y_vals <- overall_scores
          test <- calc_similarity_perm(x_vals, y_vals, method = meth, n_perm = n_perm)
          corr_val <- test$estimate
          p_val <- test$p.value
          sig <- sig_stars(p_val)
          pts <- formSymbols[as.character(combined_stats$Form)]
          plot(x_vals, y_vals,
               main = paste(score_name, "vs.", st, "\n", meth, "corr =", round(corr_val, 3), "(", sig, ")"),
               xlab = st, ylab = score_name,
               pch = pts, col = "darkgreen")
          legend("topleft", legend = paste(meth, "=", round(corr_val, 3), "(", sig, ")"), bty = "n")
          legend("bottomright", legend = names(formSymbols), pch = formSymbols, bty = "n")
        }
        base_caption(paste(setup, "All shapes pooled (symbols). PPI", score_name, "against each summary statistic; each point a genotype;", perm_txt(meth)),
                     scatter_reading)
        par(mfrow = c(1, 1))
      }
      
      ## --- Correlation Matrix Table ---
      cor_mat <- matrix(NA_real_, nrow = length(stat_names), ncol = length(score_names),
                        dimnames = list(stat_names, score_names))
      pval_mat <- matrix(NA_real_, nrow = length(stat_names), ncol = length(score_names),
                         dimnames = list(stat_names, score_names))
      for (sc_idx in seq_along(score_names)) {
        sc_name <- score_names[sc_idx]
        score_df <- scores_list[[sc_name]]
        overall_scores <- c()
        for (form in genotype_forms) {
          score_subset <- as.numeric(as.character(prepare_score_vector(score_df, form)))
          overall_scores <- c(overall_scores, score_subset)
        }
        for (i in seq_along(stat_names)) {
          st <- stat_names[i]
          x_vals <- combined_stats[[st]]
          test <- calc_similarity_perm(overall_scores, x_vals, method = meth, n_perm = n_perm)
          cor_mat[i, sc_name] <- round(test$estimate, 3)
          pval_mat[i, sc_name] <- round(test$p.value, 3)
        }
      }
      
      comb_mat <- matrix(NA_character_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                         dimnames = dimnames(cor_mat))
      for (i in seq_len(nrow(cor_mat))) {
        for (j in seq_len(ncol(cor_mat))) {
          comb_mat[i, j] <- paste0(cor_mat[i, j], " (", sig_stars(pval_mat[i, j]), ")")
        }
      }
      
      grid.newpage()
      grid.text(paste("Global Correlation Matrix - Method:", meth),
                gp = gpar(fontsize = 10, fontface = "bold"), y = 0.97)
      
      tbl <- tableGrob(comb_mat, rows = rownames(cor_mat), cols = colnames(cor_mat),
                       theme = custom_theme)
      
      for (j in seq_along(score_names)) {
        sc_name <- score_names[j]
        col_vals <- cor_mat[, sc_name]
        top_indices <- order(abs(col_vals), decreasing = TRUE)[1:2]
        
        for (idx in 1:2) {
          # Get row index for this correlation
          row_idx <- top_indices[idx]
          # Find the cell explicitly by row and column
          cell_idx <- which(tbl$layout$name == "core-bg" &
                              tbl$layout$t == (row_idx + 1) &
                              tbl$layout$l == (j + 1))
          
          if (length(cell_idx) > 0) {
            tbl$grobs[[cell_idx]]$gp$fill <- if (idx == 1) "green" else "yellow"
          } else {
            message("Could not find cell for ", sc_name, " at row ", row_idx)
          }
        }
      }
      
      pushViewport(viewport(width = 0.9, height = 1))
      grid.draw(tbl)
      popViewport()
      grid_page_caption(paste(setup, "All shapes pooled.", perm_txt(meth), "Rows: summary statistics, columns: PPIs; green: strongest, yellow: second strongest absolute correlation of each PPI."),
                        table_reading)
      
      for (sc_idx in seq_along(score_names)) {
        sc_name <- score_names[sc_idx]
        col_vals <- cor_mat[, sc_name]
        top_indices <- order(abs(col_vals), decreasing = TRUE)[1:2]
        message("GLOBAL - Method: ", meth, " | Score: ", sc_name,
                " => Top correlations: ", paste(rownames(cor_mat)[top_indices], collapse = ", "))
      }
      
      ## --- Bar Plot for Average Significant Correlations ---
      sig_cor_mat <- cor_mat
      sig_cor_mat[pval_mat >= 0.05] <- NA  # set non-significant values to NA
      avg_sig_corr <- colMeans(sig_cor_mat, na.rm = TRUE)
      
      # Create a new page and plot the barplot with fixed y-axis limits from -1 to 1
      grid.newpage()
      caption_margin(4)
      bp1 <- barplot(avg_sig_corr,
                     main = paste("Average Significant Correlation for Scores (Method:", meth, ")"),
                     xlab = "Scores",
                     ylab = "Average Significant Correlation",
                     col = "skyblue",
                     ylim = c(-1, 1),   # Force y-axis range from -1 to 1
                     yaxt = "n",        # Suppress automatic y-axis
                     las = 2)
      # Add custom y-axis with ticks from -1 to 1
      axis(2, at = seq(-1, 1, by = 0.25), las = 2)
      box()
      base_caption(paste(setup, "All shapes pooled. Mean of the significant (p < 0.05)", meth, "correlations of each PPI with the nine summary statistics."),
                   "A mean near 0 can hide strong correlations of opposite sign; use the table for the individual statistics.")
      
      ## --- Bar Plot for Average Correlations (Range, Slope, Variance) ---
      subset_stats <- c("range", "slope", "variance")
      subset_cor_mat <- cor_mat[subset_stats, , drop = FALSE]
      avg_subset_corr <- colMeans(subset_cor_mat, na.rm = TRUE)
      
      grid.newpage()
      caption_margin(4)
      bp2 <- barplot(avg_subset_corr,
                     main = paste("Average Correlation (Range, Slope, Variance) for Scores (Method:", meth, ")"),
                     xlab = "Scores",
                     ylab = "Average Correlation (Range, Slope, Variance)",
                     col = "lightgreen",
                     ylim = c(-1, 1),   # Force y-axis range from -1 to 1
                     yaxt = "n",        # Suppress automatic y-axis
                     las = 2)
      # Add custom y-axis with ticks from -1 to 1
      axis(2, at = seq(-1, 1, by = 0.25), las = 2)
      box()
      base_caption(paste(setup, "All shapes pooled. Mean", meth, "correlation of each PPI with the spread and direction statistics range, slope and variance."),
                   "High values: the PPI follows how strongly the trait changes across environments rather than its level.")
    }
    dev.off()
    message("Global PDF saved as 'scores_vs_summary_stats_global.pdf'")
  }
  
  invisible(NULL)
}




# Run the analyses.
#compare_summary_stats_with_scores(
#  rn_matrix_list = rn_matrix_list,
#  scores_list = scores_list,
#  genotype_forms = genotype_forms,
#  global = FALSE,
#  measures = c("kendall")
#)
#
#compare_summary_stats_with_scores(
#  rn_matrix_list = rn_matrix_list,
#  scores_list = scores_list,
#  genotype_forms = genotype_forms,
#  global = TRUE,
#  measures = c("kendall")
#)
