# =============================================================================
# clustering_rn_summary_stats.R — cluster reaction-norm summary statistics
# =============================================================================
# WHAT IT DOES: Clusters genotypes / indices on their reaction-norm summary statistics
#   (k-means, DBSCAN, hierarchical, mclust) and renders the associated heatmaps/plots.
#   Auto-sources 03_plasticity_scores.R when scores are not already loaded
#   (guarded by the `data_loaded` flag at the top).
# REQUIRES:     03 (so HERITABILITY_5TH / CAUSAL_SNP_NUM / OUTPUT_BASE must be set).
# PRODUCES:     Clustering figures under OUTPUT_BASE.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","clustering_rn_summary_stats.R"))
# -----------------------------------------------------------------------------
# PARAMETERS: inherits 03's parameters via the auto-source; clustering knobs are set
#   in the body of the script (edit them there). Set `data_loaded <- TRUE` beforehand to
#   reuse already-computed scores instead of re-running 03.
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

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


scores_list = list(
  CV_t = CV_t,
  RN= RN,
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
  EVS = EVS
  #MVPi = MVPi, multivariate
  #Plasticity_Ratio = Plasticity_Ratio
)

linear_matrix=do.call(rbind,linear)
gaussian_matrix=do.call(rbind,gaussian)
sinusoidal_matrix=do.call(rbind,sinusoidal)
wave_matrix=do.call(rbind,wave)

rn_matrix_list = list(
  linear = linear_matrix,    
  gaussian = gaussian_matrix, 
  sinusoidal = sinusoidal_matrix,
  wave = wave_matrix
)

genotype_forms = c("linear", "gaussian", "sinusoidal", "wave")
prepare_score_vector = function(score_df, genotype_form) {
  subset(score_df, score_df[, "Type"] == genotype_form)[, "Score"]
}

calc_similarity = function(x, y, method = "pearson") {
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

compute_summary_stats = function(rn_matrix) {
  n = ncol(rn_matrix)
  if(n < 2) stop("Reaction norm matrix must have at least 2 columns.")
  
  x_vals = seq_len(n)
  slopes = apply(rn_matrix, 1, function(x) coef(lm(x ~ x_vals))[2])
  ranges = apply(rn_matrix, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  variances = apply(rn_matrix, 1, function(x) var(x, na.rm = TRUE))
  
  half = floor(n/2)
  mean_lower = apply(rn_matrix, 1, function(x) mean(x[1:half], na.rm = TRUE))
  mean_upper = apply(rn_matrix, 1, function(x) mean(x[(half+1):n], na.rm = TRUE))
  
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


compare_clustering <- function(rn_matrix_list, scores_list, genotype_forms,
                               k_scores = 2, k_summary = 2,
                               global = FALSE,
                               pdf_file = "clustering_combined.pdf",
                               return_results = FALSE) {
  pdf(pdf_file, width = 10, height = 8)
  
  ## Helper function to color table grob cells (using the "core-bg" components)
  color_table_cells <- function(mat, tbl, col1 = "blue", col2 = "yellow") {
    nr <- nrow(mat)
    nc <- ncol(mat)
    
    for(i in seq_len(nr)) {
      for(j in seq_len(nc)) {
        cell_value <- mat[i, j]
        if(!is.na(cell_value)) {
          # Use "core-bg" to target the background
          cell_index <- which(tbl$layout$t == i + 1 & tbl$layout$l == j + 1 & tbl$layout$name == "core-bg")
          
          if(length(cell_index) > 0) {
            idx <- cell_index[1]
            if(idx <= length(tbl$grobs) && !is.null(tbl$grobs[[idx]])) {
              if(cell_value == 1) {
                tbl$grobs[[idx]]$gp <- gpar(fill = col1, col = "black")
              } else if(cell_value == 2) {
                tbl$grobs[[idx]]$gp <- gpar(fill = col2, col = "black")
              }
            }
          }
        }
      }
    }
    return(tbl)
  }
  
  # Initialize objects to store clustering results
  plasticity_clusters <- list()
  summary_clusters <- list()
  
  if (!global) {
    message("Performing LOCAL (per genotype form) score clustering...")
    for (form in genotype_forms) {
      local_cluster_list <- list()
      par(mfrow = c(3,3))
      
      # Plasticity (score) clustering for each score type
      for (score_name in names(scores_list)) {
        score_df <- scores_list[[score_name]]
        # Filter rows based on the form (assumed to be in the 2nd column)
        form_rows <- which(score_df[, 2] == form)
        if (length(form_rows) < 1) {
          message(paste("Skipping", score_name, "for form", form, "- no data"))
          next
        }
        
        # Extract score values from column 1 and genotype names from column 3 (if available)
        score_vec <- as.numeric(as.character(score_df[form_rows, 1]))
        if (ncol(score_df) >= 3) {
          geno_names <- as.character(score_df[form_rows, 3])
        } else {
          geno_names <- paste0("geno_", seq_along(score_vec))
        }
        names(score_vec) <- geno_names
        
        if (length(score_vec) < 2 || all(is.na(score_vec))) {
          message(paste("Skipping", score_name, "for form", form, "- insufficient data"))
          next
        }
        score_vec <- score_vec[!is.na(score_vec)]
        
        # Perform clustering
        d <- dist(matrix(score_vec, ncol = 1))
        hc <- hclust(d, method = "ward.D2")
        clusters <- cutree(hc, k = k_scores)
        local_cluster_list[[score_name]] <- clusters
        
        # Plot the clustering
        plot(1:length(score_vec), score_vec,
             col = clusters,
             pch = 19,
             main = paste(score_name, "\nForm:", form),
             xlab = "Genotype Index", ylab = "Score")
        cluster_centers <- tapply(score_vec, clusters, mean)
        for (cl in unique(clusters)) {
          points(which(clusters == cl),
                 rep(cluster_centers[as.character(cl)], sum(clusters == cl)), 
                 col = cl, pch = 8, cex = 2)
        }
        legend("topright", legend = paste("Cluster", sort(unique(clusters))),
               col = sort(unique(clusters)), pch = 19, cex = 0.8)
      }
      
      par(mfrow = c(1,1))
      if (length(local_cluster_list) > 0) {
        all_genos <- unique(unlist(lapply(local_cluster_list, names)))
        if (length(all_genos) > 0) {
          score_mat <- matrix(NA, nrow = length(all_genos), ncol = length(local_cluster_list))
          rownames(score_mat) <- all_genos
          colnames(score_mat) <- names(local_cluster_list)
          for (m in names(local_cluster_list)) {
            vec <- local_cluster_list[[m]]
            score_mat[names(vec), m] <- vec
          }
          grid.newpage()
          grid.text(paste("Score Clustering Results - Form:", form),
                    gp = gpar(fontsize = 10, fontface = "bold"), y = 0.97)
          custom_theme <- ttheme_minimal(
            core = list(fg_params = list(cex = 0.3)),
            colhead = list(fg_params = list(cex = 0.3)),
            rowhead = list(fg_params = list(cex = 0.3)))
          tbl <- tableGrob(as.data.frame(score_mat), rows = rownames(score_mat), theme = custom_theme)
          tbl <- color_table_cells(score_mat, tbl, col1 = "blue", col2 = "yellow")
          grid.draw(tbl)
          
          plasticity_clusters[[form]] <- local_cluster_list
        } else {
          message("Score clustering: No genotype names found for form ", form)
        }
      }
      
      # Summary statistics clustering for each summary measure
      {
        rn_matrix <- rn_matrix_list[[form]]
        summary_stats <- compute_summary_stats(rn_matrix)
        summary_stats$genotype <- rownames(summary_stats)
        names_vec <- summary_stats$genotype
        summary_stats <- summary_stats[, c("genotype", setdiff(names(summary_stats), "genotype"))]
        
        summary_cluster_list <- list()
        par(mfrow = c(3,3))
        
        for (st in setdiff(names(summary_stats), "genotype")) {
          measure_vals <- summary_stats[[st]]
          names(measure_vals) <- names_vec
          if (length(measure_vals) < 2) next
          
          d <- dist(measure_vals)
          hc <- hclust(d, method = "ward.D2")
          clusters <- cutree(hc, k = k_summary)
          summary_cluster_list[[st]] <- clusters
          
          plot(measure_vals,
               col = clusters,
               pch = 19,
               main = paste(st, "\nForm:", form),
               xlab = "Genotype Index", ylab = st)
          cluster_centers <- tapply(measure_vals, clusters, mean)
          center_indices <- sapply(unique(clusters), function(cl) which(clusters == cl)[1])
          points(center_indices, cluster_centers, col = unique(clusters), pch = 8, cex = 2)
          legend("topright", legend = paste("Cluster", sort(unique(clusters))),
                 col = sort(unique(clusters)), pch = 19, cex = 0.8)
        }
        par(mfrow = c(1,1))
        
        if (length(summary_cluster_list) > 0) {
          all_genos <- names_vec
          sum_mat <- matrix(NA, nrow = length(all_genos), ncol = length(summary_cluster_list))
          rownames(sum_mat) <- all_genos
          colnames(sum_mat) <- names(summary_cluster_list)
          for (m in names(summary_cluster_list)) {
            vec <- summary_cluster_list[[m]]
            sum_mat[names(vec), m] <- vec
          }
          grid.newpage()
          tbl_sum <- tableGrob(as.data.frame(sum_mat), rows = rownames(sum_mat),
                               theme = ttheme_minimal(
                                 core = list(fg_params = list(cex = 0.3)),
                                 colhead = list(fg_params = list(cex = 0.3)),
                                 rowhead = list(fg_params = list(cex = 0.3))))
          tbl_sum <- color_table_cells(sum_mat, tbl_sum, col1 = "blue", col2 = "yellow")
          grid.draw(tbl_sum)
          
          summary_clusters[[form]] <- summary_cluster_list
        }
      }
    }
  } else {  
    # GLOBAL clustering mode
    message("Performing GLOBAL (across genotype forms) score clustering...")
    global_cluster_list <- list()
    all_genos <- c()
    par(mfrow = c(3,3))
    
    for (score_name in names(scores_list)) {
      score_df <- scores_list[[score_name]]
      global_vec <- c()
      for (form in genotype_forms) {
        form_rows <- which(score_df[, 2] == form)
        if (length(form_rows) < 1) next
        score_subset <- as.numeric(as.character(score_df[form_rows, 1]))
        if (ncol(score_df) >= 3) {
          geno_names <- as.character(score_df[form_rows, 3])
        } else {
          geno_names <- paste0("geno_", seq_along(score_subset))
        }
        names(score_subset) <- paste0(form, "_", geno_names)
        score_subset <- score_subset[!is.na(score_subset)]
        if (length(score_subset) < 1) next
        global_vec <- c(global_vec, score_subset)
      }
      if (length(global_vec) < 2) next
      all_genos <- union(all_genos, names(global_vec))
      d <- dist(global_vec)
      hc <- hclust(d, method = "ward.D2")
      clusters <- cutree(hc, k = k_scores)
      global_cluster_list[[score_name]] <- clusters
      
      plot(global_vec,
           col = clusters,
           pch = 19,
           main = paste(score_name, "\nGlobal Clustering"),
           xlab = "Genotype Index", ylab = "Score")
      cluster_centers <- tapply(global_vec, clusters, mean)
      center_indices <- sapply(unique(clusters), function(cl) which(clusters == cl)[1])
      points(center_indices, cluster_centers, col = unique(clusters), pch = 8, cex = 2)
      legend("topright", legend = paste("Cluster", sort(unique(clusters))),
             col = sort(unique(clusters)), pch = 19, cex = 0.8)
    }
    par(mfrow = c(1,1))
    
    if (length(global_cluster_list) > 0) {
      score_mat <- matrix(NA, nrow = length(all_genos), ncol = length(global_cluster_list))
      rownames(score_mat) <- all_genos
      colnames(score_mat) <- names(global_cluster_list)
      for (m in names(global_cluster_list)) {
        vec <- global_cluster_list[[m]]
        score_mat[names(vec), m] <- vec
      }
      
      max_rows_per_page <- 30
      num_pages <- ceiling(nrow(score_mat) / max_rows_per_page)
      
      for(page in seq_len(num_pages)) {
        grid.newpage()
        start_row <- (page-1)*max_rows_per_page + 1
        end_row <- min(page*max_rows_per_page, nrow(score_mat))
        score_mat_page <- score_mat[start_row:end_row, , drop = FALSE]
        
        grid.text(paste("Global Score Clustering (Page", page, "of", num_pages, ")"),
                  gp = gpar(fontsize = 8, fontface = "bold"), y = 0.98)
        
        tbl_global <- tableGrob(as.data.frame(score_mat_page), 
                                rows = rownames(score_mat_page),
                                theme = ttheme_minimal(
                                  core = list(fg_params = list(cex = 0.35)),
                                  colhead = list(fg_params = list(cex = 0.35)),
                                  rowhead = list(fg_params = list(cex = 0.35))))
        tbl_global <- color_table_cells(score_mat_page, tbl_global, col1 = "blue", col2 = "yellow")
        grid.draw(tbl_global)
      }
      plasticity_clusters <- global_cluster_list  # for global mode, save plasticity clusters
    }
    
    global_summary <- do.call(rbind, lapply(genotype_forms, function(form) {
      rn_matrix <- rn_matrix_list[[form]]
      stats <- compute_summary_stats(rn_matrix)
      stats$genotype <- paste(form, rownames(stats), sep = "_")
      stats
    }))
    global_summary <- global_summary[, c("genotype", setdiff(names(global_summary), "genotype"))]
    names_vec <- global_summary$genotype
    
    global_summary_cluster_list <- list()
    
    par(mfrow = c(3,3))
    
    for (st in setdiff(names(global_summary), "genotype")) {
      measure_vals <- global_summary[[st]]
      names(measure_vals) <- names_vec
      if (length(measure_vals) < 2) next
      
      d <- dist(measure_vals)
      hc <- hclust(d, method = "ward.D2")
      clusters <- cutree(hc, k = k_summary)
      global_summary_cluster_list[[st]] <- clusters
      
      plot(measure_vals,
           col = clusters,
           pch = 19,
           main = paste(st, "\nGlobal Summary Clustering"),
           xlab = "Genotype Index", ylab = st)
      cluster_centers <- tapply(measure_vals, clusters, mean)
      center_indices <- sapply(unique(clusters), function(cl) which(clusters == cl)[1])
      points(center_indices, cluster_centers, col = unique(clusters), pch = 8, cex = 2)
      legend("topright", legend = paste("Cluster", sort(unique(clusters))),
             col = sort(unique(clusters)), pch = 19, cex = 0.8)
    }
    par(mfrow = c(1,1))
    
    if (length(global_summary_cluster_list) > 0) {
      sum_mat <- matrix(NA, nrow = length(names_vec), ncol = length(global_summary_cluster_list))
      rownames(sum_mat) <- names_vec
      colnames(sum_mat) <- names(global_summary_cluster_list)
      for (m in names(global_summary_cluster_list)) {
        vec <- global_summary_cluster_list[[m]]
        sum_mat[names(vec), m] <- vec
      }
      
      max_rows_per_page <- 30
      num_pages <- ceiling(nrow(sum_mat) / max_rows_per_page)
      
      for(page in seq_len(num_pages)) {
        grid.newpage()
        start_row <- (page-1)*max_rows_per_page + 1
        end_row <- min(page*max_rows_per_page, nrow(sum_mat))
        sum_mat_page <- sum_mat[start_row:end_row, , drop = FALSE]
        
        grid.text(paste("Global Summary Clustering (Page", page, "of", num_pages, ")"),
                  gp = gpar(fontsize = 8, fontface = "bold"), y = 0.98)
        
        tbl_global_sum <- tableGrob(as.data.frame(sum_mat_page), 
                                    rows = rownames(sum_mat_page),
                                    theme = ttheme_minimal(
                                      core = list(fg_params = list(cex = 0.35)),
                                      colhead = list(fg_params = list(cex = 0.35)),
                                      rowhead = list(fg_params = list(cex = 0.35))))
        tbl_global_sum <- color_table_cells(sum_mat_page, tbl_global_sum, col1 = "blue", col2 = "yellow")
        grid.draw(tbl_global_sum)
      }
      summary_clusters <- global_summary_cluster_list  # for global mode, save summary clusters
    }
  }
  
  dev.off()
  message("PDF saved as ", pdf_file)
  
  if (return_results) {
    return(list(plasticity = plasticity_clusters, summary = summary_clusters))
  }
  
  invisible(NULL)
}

##############################################
# Function to Compute ARI Between Clusterings  #
##############################################
compute_ARI_matrix <- function(plasticity_clusters, summary_clusters) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("The 'mclust' package is required. Please install it using install.packages('mclust').")
  }
  
  ari_mat <- matrix(NA, nrow = length(plasticity_clusters), ncol = length(summary_clusters))
  rownames(ari_mat) <- names(plasticity_clusters)
  colnames(ari_mat) <- names(summary_clusters)
  
  for (i in seq_along(plasticity_clusters)) {
    for (j in seq_along(summary_clusters)) {
      cl1 <- plasticity_clusters[[i]]
      cl2 <- summary_clusters[[j]]
      common_genos <- intersect(names(cl1), names(cl2))
      if (length(common_genos) < 2) {
        ari_mat[i, j] <- NA
      } else {
        ari_mat[i, j] <- mclust::adjustedRandIndex(cl1[common_genos], cl2[common_genos])
      }
    }
  }
  return(ari_mat)
}


# Run the clustering function (here in local mode) and return the clustering results
results <- compare_clustering(rn_matrix_list, scores_list, genotype_forms,
                              k_scores = 2, k_summary = 2,
                              global = FALSE,
                              pdf_file = "clustering_combined_local.pdf",
                              return_results = TRUE)

# For a specific genotype form (e.g., "form1"), extract the clustering results:
plasticity_clusters <- results$plasticity[["form1"]]
summary_clusters <- results$summary[["form1"]]
print(names(results$plasticity))

# Compute the ARI matrix comparing plasticity clusterings vs. summary clusterings for "form1"
ari_matrix <- compute_ARI_matrix(plasticity_clusters, summary_clusters)
print(ari_matrix)



global=compare_clustering(
  rn_matrix_list = rn_matrix_list,
  scores_list = scores_list,
  genotype_forms = genotype_forms,
  k_scores = 2,
  k_summary = 2,
  global = TRUE,
  pdf_file = "clustering_combined_global.pdf"
)

print(global)
