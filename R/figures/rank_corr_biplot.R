# Script: comprehensive_plasticity_analysis_corrected.R
# Purpose: Comprehensive analysis of plasticity scores with consistent layout
# All correlations use the same clustering and visualization

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(psych)
  library(reshape2)
  library(tidygraph)
  library(ggraph)
  library(grid)
  library(scales)
  library(tibble)
})




is_maize <- exists("form_ranges") && length(form_ranges) == 3


# Which maize trait to analyse. NULL/"all" = pool all 732 rows across biomass +
# leafArea + WUE (mirrors the synthetic run, which pools all four forms).
# Set to one of names(form_ranges) ("biomass", "leafArea", "WUE") to subset.
if (is_maize && !exists("maize_trait")) {
  maize_trait <- "biomass"
}
maize_pool_all <- !is_maize || is.null(maize_trait) || identical(maize_trait, "all")

if (!exists("output_dir")) {
  output_dir <- if (is_maize) {
    suffix <- if (maize_pool_all) "" else paste0("_", maize_trait)
    file.path("~/CRC_1644_Z2_GWAS_simple/plasticity_comprehensive_analysis_maize")
  } else {
    "~/CRC_1644_Z2_GWAS_simple/plasticity_comprehensive_analysis"
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

regression_data_path <- "~/CRC_1644_Z2_GWAS_simple/R-files/regression_summary_stats/regression_data_full_interval_10_indices_1_11_21_31_41_50.csv"
maize_scores_path    <- "~/CRC_1644_Z2_GWAS_simple/maize_scores_list.rds"

tau_threshold <- 0.8
topk_values <- c(3, 5, 10, 15, 20)
AGREE_THR_FULL <- 1.00

sumstats_ps_cols <- 2:29
sumstats_stat_from <- 30
# ======================================================

# ===================== UTILITIES ======================
theme_pub <- function(base_size = 10, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(plot.title.position = "plot",
          plot.caption.position = "plot",
          axis.title = element_text(),
          axis.text = element_text(),
          legend.position = "right",
          legend.title = element_text(),
          panel.grid.major = element_line(linewidth = 0.25),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
}

scale_fill_div_pub <- function(lim = c(-1,1)) {
  scale_fill_gradient2(limits = lim, oob = squish, midpoint = 0,
                       low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426")
}

scale_color_disc_pub <- function() scale_color_manual(values = scales::hue_pal(l = 55)(10))
scale_fill_disc_pub <- function() scale_fill_manual(values = scales::hue_pal(l = 55)(10))

load_regression_data <- function(data_path) {
  if (!file.exists(data_path)) {
    stop("Regression data file not found: ", data_path)
  }
  
  df <- read_csv(data_path, show_col_types = FALSE)
  
  plasticity_scores <- df[, sumstats_ps_cols] %>%
    dplyr::select(where(is.numeric)) %>%       # <--- Added dplyr::
    as.matrix()
  
  summary_stats <- df[, sumstats_stat_from:ncol(df)] %>%
    dplyr::select(where(is.numeric)) %>%       # <--- Added dplyr::
    as.matrix()
  
  if ("...1" %in% colnames(df)) {
    rownames(plasticity_scores) <- df$`...1`
    rownames(summary_stats) <- df$`...1`
  } else {
    rownames(plasticity_scores) <- paste0("Genotype_", 1:nrow(plasticity_scores))
    rownames(summary_stats) <- paste0("Genotype_", 1:nrow(summary_stats))
  }
  
  return(list(
    plasticity_scores = plasticity_scores,
    summary_stats = summary_stats,
    genotype_names = rownames(plasticity_scores)
  ))
}

# CORRECTED: Proper correlation computations
compute_correlations <- function(plasticity_scores) {
  # Kendall and Spearman on raw scores (they handle ranking internally)
  ct_tau <- corr.test(plasticity_scores, method = "kendall", adjust = "none")
  ct_spr <- corr.test(plasticity_scores, method = "spearman", adjust = "none")
  
  # Pearson on absolute values
  ct_pea_abs <- corr.test(plasticity_scores, method = "pearson", adjust = "none")
  
  # For top-k analysis: ranks (higher score = better rank)
  rank_matrix <- apply(plasticity_scores, 2, function(x) rank(-x, ties.method = "average"))
  
  return(list(
    kendall = ct_tau$r,
    spearman = ct_spr$r,
    pearson_absolute = ct_pea_abs$r,
    p_kendall = ct_tau$p,
    p_spearman = ct_spr$p,
    p_pearson_absolute = ct_pea_abs$p,
    ranks = rank_matrix
  ))
}

compute_method_stat_correlations <- function(plasticity_scores, summary_stats) {
  combined_data <- cbind(plasticity_scores, summary_stats)
  ct_ms <- corr.test(combined_data, method = "spearman", adjust = "none")
  
  r_ms <- ct_ms$r[colnames(plasticity_scores), colnames(summary_stats)]
  p_ms <- ct_ms$p[colnames(plasticity_scores), colnames(summary_stats)]
  
  cor_df <- reshape2::melt(r_ms, varnames = c("PS", "Stat"), value.name = "r") %>%
    left_join(
      reshape2::melt(p_ms, varnames = c("PS", "Stat"), value.name = "p"),
      by = c("PS", "Stat")
    ) %>%
    mutate(star = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ))
  
  return(list(
    cor_df = cor_df,
    r_ms = r_ms
  ))
}
# ======================================================

# ================= CORE CALCULATIONS ==================
topk_set_idx <- function(ranks_df, method, k) {
  v <- ranks_df[[method]]
  ord <- order(v, na.last = NA)
  head(ranks_df$gid[ord], min(k, sum(!is.na(v))))
}

topk_overlap_tbl <- function(ranks_df, ks) {
  idxs <- setdiff(colnames(ranks_df), "gid")
  if (length(idxs) < 2) return(tibble())
  
  pairs <- combn(idxs, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pr) {
    a <- pr[1]; b <- pr[2]
    A_all <- ranks_df$gid[order(ranks_df[[a]], na.last = NA)]
    B_all <- ranks_df$gid[order(ranks_df[[b]], na.last = NA)]
    
    tibble(score_i = a, score_j = b, k = ks) |>
      mutate(overlap = purrr::map_int(k, ~length(intersect(
        A_all[seq_len(min(.x, length(A_all)))],
        B_all[seq_len(min(.x, length(B_all)))]
      ))),
      overlap_frac = overlap / k)
  })
}

topk_overlap_matrix <- function(ranks_df, methods, k) {
  m <- length(methods)
  M <- matrix(NA_real_, m, m, dimnames = list(methods, methods))
  sets <- lapply(methods, function(md) topk_set_idx(ranks_df, md, k))
  
  for (i in seq_len(m)) {
    for (j in i:m) {
      ov <- length(intersect(sets[[i]], sets[[j]])) / k
      M[i, j] <- M[j, i] <- ov
    }
  }
  M
}
# ======================================================

# ================= PLOTTING FUNCTIONS =================
# Use EXACT same functions from your original script for consistency
plot_tau_heatmap <- function(Tm, tag) {
  m <- Tm; diag(m) <- 1; m[is.na(m)] <- 0
  d <- as.dist(1 - m); hc <- hclust(d, method = "average")
  ord <- rownames(m)[hc$order]
  df <- as.data.frame(m) |>
    tibble::rownames_to_column("score_i") |>
    tidyr::pivot_longer(-score_i, names_to = "score_j", values_to = "tau") |>
    mutate(score_i = factor(score_i, levels = ord),
           score_j = factor(score_j, levels = ord))
  ggplot(df, aes(score_j, score_i, fill = tau)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_div_pub(lim = c(-1,1)) +
    labs(x = NULL, y = NULL, fill = "τ") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
}

plot_agreement_network <- function(Tm, tau_threshold, tag) {
  A <- ifelse(Tm >= tau_threshold, 1, 0); diag(A) <- 0
  G <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  if (gorder(G) == 0) return(ggplot() + theme_void())
  comps <- components(G)$membership
  deg   <- degree(G)
  g <- as_tbl_graph(G) |> activate(nodes) |> mutate(component = factor(comps), deg = deg)
  
  ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.4) +
    geom_node_point(aes(fill = component, size = deg), shape = 21, color = "black", stroke = 0.2) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 50) +
    scale_size_continuous(range = c(3,10), breaks = scales::pretty_breaks(3)) +
    scale_fill_disc_pub() +
    guides(size = guide_legend(title = "Degree"), fill = guide_legend(title = "Component")) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(legend.key.height = grid::unit(0.5, "cm"), legend.key.width = grid::unit(0.5, "cm"))
}

plot_method_stat_heatmap_stars <- function(cor_df, method_order = NULL, tag = "") {
  df <- cor_df
  # Order the Y-axis (PS) based on the supplied clustered order
  if (!is.null(method_order)) df <- df %>% mutate(PS = factor(PS, levels = method_order))
  ggplot(df, aes(x = Stat, y = PS, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_div_pub(lim = c(-1,1)) +
    geom_text(aes(label = sprintf("%.2f%s", r, star)), size = 3) +
    labs(x = "Summary Statistic", y = "Plasticity Score / Method", fill = "τ") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
}

# Use the SAME heatmap function for ALL correlations (with clustering)
plot_consistent_corr_heatmap <- function(M, tag = "", title = "", fill_label = "corr") {
  # Apply the same clustering as plot_tau_heatmap
  m <- M; diag(m) <- 1; m[is.na(m)] <- 0
  d <- as.dist(1 - m); hc <- hclust(d, method = "average")
  ord <- rownames(m)[hc$order]
  
  df <- as.data.frame(m) |>
    tibble::rownames_to_column("method_i") |>
    tidyr::pivot_longer(-method_i, names_to = "method_j", values_to = "val") |>
    mutate(method_i = factor(method_i, levels = ord),
           method_j = factor(method_j, levels = ord))
  
  ggplot(df, aes(method_j, method_i, fill = val)) +
    geom_tile(color = "white") +
    coord_fixed() +
    scale_fill_div_pub(lim = c(-1,1)) +
    labs(x = NULL, y = NULL, fill = fill_label) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

plot_topk_overlap_heatmap_matrix <- function(Om, tag, k) {
  m <- Om; diag(m) <- 1; m[is.na(m)] <- 0
  d <- as.dist(1 - m); hc <- hclust(d, method = "average"); ord <- rownames(m)[hc$order]
  df <- as.data.frame(m) |>
    tibble::rownames_to_column("method_i") |>
    tidyr::pivot_longer(-method_i, names_to = "method_j", values_to = "overlap") |>
    mutate(method_i = factor(method_i, levels = ord),
           method_j = factor(method_j, levels = ord))
  ggplot(df, aes(method_j, method_i, fill = overlap)) +
    geom_tile() + coord_fixed() +
    scale_fill_viridis_c(limits = c(0,1), labels = percent) +
    labs(x = NULL, y = NULL, fill = "overlap") +
    theme_pub() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

# PCA Biplot using ranks
plot_genotype_pca_biplot <- function(plasticity_ranks, summary_stats, Tm, tau_threshold, tag) {
  complete_cases <- complete.cases(plasticity_ranks, summary_stats)
  ps_clean <- plasticity_ranks[complete_cases, , drop = FALSE]
  ss_clean <- summary_stats[complete_cases, , drop = FALSE]
  
  if (nrow(ps_clean) < 3) {
    return(ggplot() + theme_void() + 
             labs(title = paste0("Not enough complete genotypes for PCA • ", tag)))
  }
  
  ps_clean <- ps_clean[, apply(ps_clean, 2, function(x) length(unique(x)) > 1), drop = FALSE]
  ss_clean <- ss_clean[, apply(ss_clean, 2, function(x) length(unique(x)) > 1), drop = FALSE]
  
  if (ncol(ps_clean) < 2) {
    return(ggplot() + theme_void() + 
             labs(title = paste0("Not enough variance in plasticity ranks • ", tag)))
  }
  
  pca_result <- prcomp(ps_clean, center = TRUE, scale. = TRUE)
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  genotype_scores <- as.data.frame(pca_result$x[, 1:2])
  colnames(genotype_scores) <- c("PC1", "PC2")
  genotype_scores$genotype <- rownames(ps_clean)
  
  method_loadings <- as.data.frame(pca_result$rotation[, 1:2])
  colnames(method_loadings) <- c("PC1", "PC2")
  method_loadings$method <- colnames(ps_clean)
  method_loadings$type <- "Plasticity Method"
  
  if (ncol(ss_clean) > 0) {
    stat_correlations <- cor(ss_clean, pca_result$x[, 1:2], use = "pairwise.complete.obs")
    stat_loadings <- as.data.frame(stat_correlations)
    colnames(stat_loadings) <- c("PC1", "PC2")
    stat_loadings$method <- rownames(stat_correlations)
    stat_loadings$type <- "Summary Statistic"
  } else {
    stat_loadings <- data.frame(PC1 = numeric(), PC2 = numeric(), 
                                method = character(), type = character())
  }
  
  all_loadings <- bind_rows(method_loadings, stat_loadings)
  
  plasticity_methods <- all_loadings %>% filter(type == "Plasticity Method")
  common_methods <- intersect(rownames(Tm), plasticity_methods$method)
  
  if (length(common_methods) > 0) {
    A <- ifelse(Tm[common_methods, common_methods] >= tau_threshold, 1, 0)
    diag(A) <- 0
    Gtau <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    
    if (gorder(Gtau) > 0) {
      comps <- components(Gtau)$membership
    } else {
      comps <- setNames(rep(1, length(common_methods)), common_methods)
    }
    
    all_loadings <- all_loadings %>%
      left_join(
        data.frame(method = names(comps), component = factor(unname(comps))),
        by = "method"
      )
  } else {
    all_loadings$component <- factor(1)
  }
  
  genotype_range <- max(abs(range(genotype_scores$PC1, genotype_scores$PC2)))
  
  plasticity_loadings <- all_loadings %>% filter(type == "Plasticity Method")
  if (nrow(plasticity_loadings) > 0) {
    plasticity_range <- max(abs(range(plasticity_loadings$PC1, plasticity_loadings$PC2)))
    if (plasticity_range > 0) {
      scale_factor <- 0.7 * genotype_range / plasticity_range
      all_loadings <- all_loadings %>%
        mutate(PC1 = ifelse(type == "Plasticity Method", PC1 * scale_factor, PC1),
               PC2 = ifelse(type == "Plasticity Method", PC2 * scale_factor, PC2))
    }
  }
  
  stat_loadings <- all_loadings %>% filter(type == "Summary Statistic")
  if (nrow(stat_loadings) > 0) {
    stat_range <- max(abs(range(stat_loadings$PC1, stat_loadings$PC2)))
    if (stat_range > 0) {
      scale_factor <- 0.7 * genotype_range / stat_range
      all_loadings <- all_loadings %>%
        mutate(PC1 = ifelse(type == "Summary Statistic", PC1 * scale_factor, PC1),
               PC2 = ifelse(type == "Summary Statistic", PC2 * scale_factor, PC2))
    }
  }
  
  p <- ggplot() +
    geom_point(data = genotype_scores, 
               aes(x = PC1, y = PC2), 
               alpha = 0.5, size = 2, color = "gray60") +
    
    geom_segment(data = all_loadings %>% filter(type == "Summary Statistic"),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "inches"), type = "closed"),
                 color = "red", linewidth = 1, alpha = 0.8) +
    geom_text_repel(data = all_loadings %>% filter(type == "Summary Statistic"),
                    aes(x = PC1, y = PC2, label = method),
                    color = "red", fontface = "bold", size = 4,
                    max.overlaps = 20, box.padding = 0.5) +
    
    geom_segment(data = all_loadings %>% filter(type == "Plasticity Method"),
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, color = component),
                 arrow = arrow(length = unit(0.15, "inches")),
                 linewidth = 0.8, alpha = 0.8) +
    geom_text_repel(data = all_loadings %>% filter(type == "Plasticity Method"),
                    aes(x = PC1, y = PC2, label = method, color = component),
                    size = 3.5, max.overlaps = 20, show.legend = FALSE,
                    box.padding = 0.5) +
    
    scale_color_manual(values = scales::hue_pal(l = 55)(10), name = "τ-component") +
    labs(
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme_pub() +
    coord_equal()
  
  return(p)
}
# ======================================================

# ================= MAIN ANALYSIS ======================
main <- function() {
  cat("Starting comprehensive plasticity analysis...\n")
  
  # Load data
  if (is_maize) {
    trait_label <- if (maize_pool_all) "all (pooled)" else maize_trait
    cat("Loading maize plasticity scores from RDS (trait: ", trait_label, ")\n", sep = "")
    scores_list <- readRDS(maize_scores_path)
    scores_arr  <- array(unlist(scores_list), dim = c(732, 3, 28))
    plasticity_scores <- as.numeric(scores_arr[, 1, ])
    dim(plasticity_scores)      <- c(732, 28)
    colnames(plasticity_scores) <- names(scores_list)
    rownames(plasticity_scores) <- paste0("Genotype_", seq_len(732))

    # Optionally subset to a single trait (rows are biomass | leafArea | WUE stacked).
    # Default behaviour: pool all 732 rows, mirroring the synthetic run which pools all forms.
    if (!maize_pool_all) {
      if (!maize_trait %in% names(form_ranges)) {
        stop("maize_trait '", maize_trait, "' not found in form_ranges (",
             paste(names(form_ranges), collapse = ", "), ")")
      }
      plasticity_scores <- plasticity_scores[form_ranges[[maize_trait]], , drop = FALSE]
    }

    # No summary statistics available for the maize dataset.
    summary_stats <- matrix(numeric(0), nrow = nrow(plasticity_scores), ncol = 0,
                            dimnames = list(rownames(plasticity_scores), NULL))
  } else {
    cat("Loading regression data...\n")
    data <- load_regression_data(regression_data_path)
    plasticity_scores <- data$plasticity_scores
    summary_stats <- data$summary_stats
  }

  cat("Data loaded:\n")
  cat("  Plasticity scores:", nrow(plasticity_scores), "genotypes ×", ncol(plasticity_scores), "methods\n")
  cat("  Summary stats:", nrow(summary_stats), "genotypes ×", ncol(summary_stats), "statistics\n")

  # Compute correlations
  cat("Computing correlations...\n")
  cor_results <- compute_correlations(plasticity_scores)
  Tm <- cor_results$kendall  # Use Kendall for network and coloring

  method_stat_cor <- NULL
  if (ncol(summary_stats) > 0) {
    cat("Computing method-stat correlations...\n")
    method_stat_cor <- compute_method_stat_correlations(plasticity_scores, summary_stats)
  } else {
    cat("Skipping method-stat correlations (no summary statistics available).\n")
  }
  
  # Prepare ranks data for top-k analysis
  ranks_df <- as.data.frame(cor_results$ranks) %>%
    mutate(gid = row_number()) %>%
    relocate(gid, .before = 1)
  
  # Create output PDF
  pdf_path <- file.path(output_dir, "comprehensive_plasticity_analysis.pdf")
  pdf(pdf_path, width = 6.3, height = 7)
  
  cat("Creating analysis plots...\n")
  
  # --- Calculate the clustered order from the PS-PS Kendall correlation matrix (Tm) ---
  # This order will be reused for all heatmaps for consistency.
  m_ord <- Tm
  diag(m_ord) <- 1
  m_ord[is.na(m_ord)] <- 0
  d_ord <- as.dist(1 - m_ord)
  hc_ord <- hclust(d_ord, method = "average")
  clustering_order <- rownames(m_ord)[hc_ord$order]
  # -----------------------------------------------------------------------------------
  
  # 1. Kendall τ Heatmap (with clustering)
  cat("  Creating Kendall τ heatmap...\n")
  p1 <- plot_tau_heatmap(cor_results$kendall, "Plasticity Analysis")
  print(p1)
  
  # 2. All Correlation Heatmaps (ALL with same clustering)
  cat("  Creating correlation heatmaps...\n")
  # These functions already use the internal clustering logic, but we include them here for context
  p2 <- plot_consistent_corr_heatmap(cor_results$kendall, "Plasticity Analysis",
                                     "Kendall correlation (τ) - Rank Agreement", "τ")
  print(p2)
  
  p3 <- plot_consistent_corr_heatmap(cor_results$spearman, "Plasticity Analysis",
                                     "Spearman correlation (ρ) - Rank Agreement", "ρ")
  print(p3)
  
  p4 <- plot_consistent_corr_heatmap(cor_results$pearson_absolute, "Plasticity Analysis",
                                     "Pearson correlation (r) - Absolute Values", "r")
  print(p4)
  
  # 3. Method × Stat Correlation Heatmap (FIXED: now uses the clustered order)
  if (!is.null(method_stat_cor)) {
    cat("  Creating method-stat correlation heatmap...\n")
    p5 <- plot_method_stat_heatmap_stars(method_stat_cor$cor_df,
                                         method_order = clustering_order, # <--- FIX: Passing the clustered order
                                         tag = "Plasticity Analysis")
    print(p5)
  }
  
  # 4. Agreement Network
  cat("  Creating agreement network...\n")
  p6 <- plot_agreement_network(Tm, tau_threshold, "Plasticity Analysis")
  print(p6)
  
  # 5. Top-k Overlap Heatmaps
  cat("  Creating top-k overlap heatmaps...\n")
  for (k in topk_values) {
    Om <- topk_overlap_matrix(ranks_df, colnames(plasticity_scores), k)
    p_topk <- plot_topk_overlap_heatmap_matrix(Om, paste0("Plasticity Analysis (k=", k, ")"), k)
    print(p_topk)
  }
  
  # 6. PCA Biplot (using ranks)
  cat("  Creating PCA biplot...\n")
  p7 <- plot_genotype_pca_biplot(cor_results$ranks, summary_stats, Tm, tau_threshold, 
                                 "Plasticity Analysis")
  print(p7)
  
  dev.off()
  
  # Save data files
  cat("Saving data files...\n")
  
  kendall_df <- as.data.frame(cor_results$kendall) %>% rownames_to_column("method")
  spearman_df <- as.data.frame(cor_results$spearman) %>% rownames_to_column("method")
  pearson_abs_df <- as.data.frame(cor_results$pearson_absolute) %>% rownames_to_column("method")
  ranks_df_out <- as.data.frame(cor_results$ranks) %>% rownames_to_column("genotype")
  
  write_csv(kendall_df, file.path(output_dir, "kendall_correlation_matrix.csv"))
  write_csv(spearman_df, file.path(output_dir, "spearman_correlation_matrix.csv"))
  write_csv(pearson_abs_df, file.path(output_dir, "pearson_absolute_correlation_matrix.csv"))
  write_csv(ranks_df_out, file.path(output_dir, "plasticity_score_ranks.csv"))
  if (!is.null(method_stat_cor)) {
    write_csv(method_stat_cor$cor_df, file.path(output_dir, "method_stat_correlations.csv"))
  }
  
  cat("\nComprehensive analysis complete!\n")
  cat("Output saved to:", pdf_path, "\n")
  cat("Data files saved to:", output_dir, "\n")
}

# Run the analysis

main()