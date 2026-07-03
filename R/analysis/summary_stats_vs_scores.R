if (!exists("data_loaded") || !data_loaded) {
  source(here::here("R", "03_plasticity_scores.R"))
  data_loaded = TRUE
}

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
scores_list = list(
  CV_t = CV_t,
  RN = RN,
  RNN = RNN,
  D_slope = D_slope,
  RC = RC,
  #CVm = CVm,comparative
  #CVmd = CVmd,comparative
  gPi = gPi,
  PPF = PPF,
  PPi = PPi,
  PImd = PImd,
  PILSM = PILSM,
  RTR = RTR,
  PIR = PIR,
  RDPI = RDPI,
  #RDPI_mean = RDPI_mean, comparative 
  ESPI = ESPI,
  ESPIID = ESPIID,
  PSI = PSI,
  RPI = RPI,
  PQ = PQ,
  PR = PR,
  NRW = NRW,
  ESP = ESP,
  #PD = PD, control stress
  #FPI = FPI,control stress
  #TPS = TSP,control stress
  CEV = CEV,
  PFI = PFI,
  APC = APC,
  SI = SI,
  RSI = RSI,
  EVS = EVS,
  #MVPi = MVPi, multivariate
  #Plasticity_Ratio = Plasticity_Ratio
  FW = FW
)

# Combine matrices for each genotype form.
linear_matrix   = do.call(rbind, linear)
gaussian_matrix = do.call(rbind, gaussian)
sinusoidal_matrix = do.call(rbind, sinusoidal)
wave_matrix     = do.call(rbind, wave)

rn_matrix_list <- list(
  linear = linear_matrix,    
  gaussian = gaussian_matrix, 
  sinusoidal = sinusoidal_matrix,
  wave = wave_matrix
)

genotype_forms = c("linear", "gaussian", "sinusoidal", "wave")



prepare_score_vector <- function(score_df, genotype_form) {
  # Assumes score_df has columns "Score" and "Type"
  subset(score_df, score_df[, "Type"] == genotype_form)[, "Score"]
}

# Original similarity function using various methods.
calc_similarity <- function(x, y, method = "pearson") {
  if (method == "pearson") {
    return(cor(x, y, method = "pearson", use = "complete.obs"))
  } else if (method == "spearman") {
    return(cor(x, y, method = "spearman", use = "complete.obs"))
  } else if (method == "kendall") {
    return(cor(x, y, method = "kendall", use = "complete.obs"))
  } else if (method == "dcor") {
    return(dcor(x, y))
  } else if (method == "cosine") {
    return(sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
  } else {
    stop("Unknown method: ", method)
  }
}

calc_similarity_perm <- function(x, y, method = "pearson", n_perm = 999) {
  # Convert inputs to numeric (if not already)
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Determine indices with complete cases
  complete_idx <- complete.cases(x, y)
  if (sum(complete_idx) < 2) {
    return(list(estimate = NA, p.value = NA))
  }
  
  # Use only complete pairs
  x_complete <- x[complete_idx]
  y_complete <- y[complete_idx]
  
  # Compute the observed similarity measure
  obs <- calc_similarity(x_complete, y_complete, method)
  
  # Perform permutation test by shuffling y_complete
  perm_stats <- replicate(n_perm, {
    permuted_y <- sample(y_complete, length(y_complete), replace = FALSE)
    calc_similarity(x_complete, permuted_y, method)
  })
  
  # Two-sided test: count how many permuted similarities are as extreme as observed.
  p_val <- (sum(abs(perm_stats) >= abs(obs)) + 1) / (n_perm + 1)
  
  return(list(estimate = obs, p.value = p_val))
}

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
                                              measures = c("kendall"), n_perm = 999) {
  # Define the nine summary statistics.
  stat_names <- c("min", "max", "mean", "median", "slope", "range", "variance", "mean_lower", "mean_upper")
  score_names <- names(scores_list)
  
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
    pdf("scores_vs_summary_stats_local.pdf", width = 10, height = 8)
    
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
          
          par(mfrow = c(3, 3))
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
    
    pdf("scores_vs_summary_stats_global.pdf", width = 12, height = 8)
    
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
        par(mfrow = c(3, 3))
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
      
      ## --- Bar Plot for Average Correlations (Range, Slope, Variance) ---
      subset_stats <- c("range", "slope", "variance")
      subset_cor_mat <- cor_mat[subset_stats, , drop = FALSE]
      avg_subset_corr <- colMeans(subset_cor_mat, na.rm = TRUE)
      
      grid.newpage()
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
