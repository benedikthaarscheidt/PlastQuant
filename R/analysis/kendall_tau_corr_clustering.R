# =============================================================================
# kendall_tau_corr_clustering.R — cluster indices by rank agreement
# =============================================================================
# WHAT IT DOES: Reads the rank-agreement ("rank flops") tables and clusters the
#   plasticity indices by their Kendall-tau agreement across top-k rankings.
# REQUIRES:     Rank-flop CSVs under output/rank_flops/ (produced upstream).
# PRODUCES:     Clustering summaries/plots of index agreement.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","kendall_tau_corr_clustering.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   in_dir          path, default output/rank_flops   input rank-agreement files          [COMMON]
#   tau_threshold   numeric, default 0.8              Kendall-tau agreement threshold      [COMMON]
#   topk_values     integer vector, default c(3,5,10,15,20)  top-k cutoffs to evaluate     [COMMON]
#   AGREE_THR_FULL  numeric, default 1.00             full-agreement threshold
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(igraph)
  library(ggplot2)
  library(scales)
  library(tidygraph)
  library(ggraph)
  library(ggrepel)
  library(grid)
  library(psych)        # corr.test()
  library(reshape2)     # melt()
  library(viridis)      # viridis colour scales
})

# ======================= CONFIG =======================
in_dir  <- here::here("output", "rank_flops")
out_dir <- file.path(in_dir, "derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tau_threshold <- 0.8
topk_values  <- c(3,5,10,15,20)
AGREE_THR_FULL <- 1.00

# Summary-stats/PS realisations file used in your other figure
sumstats_path_override <-
  here::here("output", "regression_summary_stats", "regression_data_full_interval_15_indices_1_16_31_46_50.csv")
# Adjust columns if your file layout differs:
sumstats_ps_cols   <- 2:28   # plasticity index realisations (methods)
sumstats_stat_from <- 29     # first summary-stat column
# ======================================================


# ===================== UTIL & STYLES ==================
find_files <- function(dir, pat) list.files(dir, pattern = pat, full.names = TRUE)
build_tag  <- function(path, prefix) { b <- basename(path); x <- sub(paste0("^", prefix, "_"), "", b); sub("\\.csv$", "", x) }

theme_pub <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_text(face = "bold"),
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
scale_fill_disc_pub  <- function() scale_fill_manual(values = scales::hue_pal(l = 55)(10))

ensure_gid <- function(df) {
  if (!"gid" %in% names(df)) df <- df %>% arrange(genotype_id) %>%
      mutate(gid = dplyr::row_number()) %>% relocate(gid, .before = 1)
  df
}
# ======================================================


# ================= CORE CALCULATIONS ==================
exact_flips_from_ranks <- function(r1, r2) {
  n <- length(r1); conc <- 0L; disc <- 0L
  if (n >= 2) {
    for (i in 1:(n-1)) {
      di <- r1[i] - r1[(i+1):n]
      dj <- r2[i] - r2[(i+1):n]
      s  <- di * dj
      conc <- conc + sum(s > 0)
      disc <- disc + sum(s < 0)
    }
  }
  pairs <- n*(n-1)/2
  tibble(concordant = conc, discordant = disc, n_pairs = pairs,
         flips_frac = ifelse(pairs>0, disc/pairs, NA_real_))
}

topk_set_idx <- function(ranks_df, method, k) {
  v <- ranks_df[[method]]
  ord <- order(v, na.last = NA)
  head(ranks_df$gid[ord], min(k, sum(!is.na(v))))
}

topk_overlap_tbl <- function(ranks_df, ks) {
  idxs <- setdiff(colnames(ranks_df), c("gid","genotype_id"))
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
# ======================================================


# ======================= PLOTS ========================
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
    labs(title = paste0("Kendall's τ • ", tag), x = NULL, y = NULL, fill = "τ") +
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
    labs(title = paste0("Agreement network (τ ≥ ", tau_threshold, ") • ", tag), x = NULL, y = NULL) +
    theme_pub() +
    theme(legend.key.height = grid::unit(0.5, "cm"), legend.key.width = grid::unit(0.5, "cm"))
}

plot_component_sizes <- function(Tm, tau_threshold, tag) {
  A <- ifelse(Tm >= tau_threshold, 1, 0); diag(A) <- 0
  G <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  if (gorder(G) == 0) return(ggplot() + theme_void())
  cs <- table(components(G)$membership)
  df <- tibble::tibble(component = factor(names(cs), levels = names(cs)), size = as.integer(cs))
  ggplot(df, aes(component, size, fill = component)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_disc_pub() +
    labs(title = paste0("Component sizes • ", tag), x = "Component", y = "Number of methods") +
    theme_pub() +
    theme(legend.position = "none")
}

plot_flips_vs_tau <- function(flips_tbl, pairs_df, tag) {
  df <- flips_tbl |>
    dplyr::inner_join(dplyr::select(pairs_df, score_i, score_j, correlation),
                      by = c("score_i","score_j")) |>
    dplyr::mutate(flips_pct = 100*flips_frac)
  ggplot(df, aes(correlation, flips_pct)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste0("Discordant pairs vs Kendall’s τ • ", tag),
         x = "Kendall’s τ", y = "Discordant pairs (%)") +
    theme_pub()
}

method_agreement_metrics <- function(Tm, tau_threshold, topk_df, flips_tbl) {
  idx <- rownames(Tm)
  deg <- rowSums((Tm >= tau_threshold) & !is.na(Tm)) - 1L
  avg_tau <- rowMeans(Tm, na.rm = TRUE)
  med_tau <- apply(Tm, 1, function(x) median(x, na.rm = TRUE))
  avg_topk <- purrr::map_dfr(idx, function(m) {
    df <- topk_df |> filter(score_i == m | score_j == m) |> group_by(k) |>
      summarize(avg_overlap = mean(overlap_frac, na.rm = TRUE), .groups = "drop")
    if (nrow(df) == 0) tibble(method = m, k = NA_integer_, avg_overlap = NA_real_)
    else mutate(df, method = m, .before = 1)
  })
  avg_topk_overall <- avg_topk |> group_by(method) |>
    summarize(avg_topk_overlap = mean(avg_overlap, na.rm = TRUE), .groups = "drop")
  flips_long <- flips_tbl |>
    mutate(method = score_i) |>
    bind_rows(flips_tbl |> transmute(method = score_j, flips_frac)) |>
    group_by(method) |>
    summarize(avg_flips_frac = mean(flips_frac, na.rm = TRUE), .groups = "drop")
  tibble(method = idx) |>
    mutate(degree_ge_thresh = deg, avg_tau = avg_tau, med_tau = med_tau) |>
    left_join(avg_topk_overall, by = "method") |>
    left_join(flips_long,      by = "method") |>
    arrange(desc(avg_tau))
}

plot_method_agreement <- function(metrics, tau_threshold, tag) {
  ggplot(metrics, aes(x = avg_tau, y = avg_topk_overlap, size = degree_ge_thresh)) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = tau_threshold, linetype = 2) +
    scale_size_continuous(name = paste0("deg(τ≥", tau_threshold, ")")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste0("Method agreement summary • ", tag),
         x = "Average Kendall’s τ vs others",
         y = "Average top-k overlap vs others") +
    theme_pub()
}

plot_method_flips <- function(metrics, tag) {
  ggplot(metrics, aes(x = avg_tau, y = 100*avg_flips_frac, size = degree_ge_thresh)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(name = "Graph degree") +
    labs(title = paste0("τ vs discordance • ", tag),
         x = "Average Kendall’s τ vs others",
         y = "Discordant pairs (%)") +
    theme_pub()
}

method_overlap_summary_from_topk <- function(topk_df, agree_thr = AGREE_THR_FULL) {
  if (nrow(topk_df) == 0) return(tibble())
  both <- dplyr::bind_rows(
    topk_df %>% transmute(k, method = score_i, other = score_j, overlap_frac),
    topk_df %>% transmute(k, method = score_j, other = score_i, overlap_frac)
  )
  both %>%
    group_by(method, k) %>%
    summarize(
      avg_overlap = mean(overlap_frac, na.rm = TRUE),
      min_overlap = suppressWarnings(min(overlap_frac, na.rm = TRUE)),
      frac_ge_thr = mean(overlap_frac >= agree_thr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(full_agree_all = (min_overlap >= agree_thr))
}

plot_method_k_heatmap <- function(mos, tag, value_col = "avg_overlap",
                                  title_suffix = "Average top-k overlap vs others (global)") {
  if (nrow(mos) == 0) return(ggplot() + theme_void())
  
  ord <- mos %>%
    dplyr::group_by(method) %>%
    dplyr::summarize(s = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(s)) %>% dplyr::pull(method)
  
  df <- mos %>%
    dplyr::transmute(method, k, value = .data[[value_col]], full_agree_all) %>%
    dplyr::mutate(method = factor(method, levels = ord))
  
  ggplot(df, aes(x = factor(k), y = method, fill = value)) +
    geom_tile() +
    geom_point(
      data = df[df$full_agree_all %in% TRUE, , drop = FALSE],
      aes(x = factor(k), y = method),
      shape = 21, fill = "white", size = 2, stroke = 0.3, inherit.aes = FALSE
    ) +
    scale_fill_viridis_c(limits = c(0,1), labels = scales::percent, name = value_col) +
    labs(title = paste0(title_suffix, " • ", tag), x = "k", y = "Method") +
    theme_pub()
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
    labs(title = paste0("Pairwise top-", k, " overlap (methods × methods) • ", tag),
         x = NULL, y = NULL, fill = "overlap") +
    theme_pub() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

threshold_sweep_stats <- function(Tm, tau_seq = seq(0.3, 0.9, by = 0.05)) {
  bind_rows(lapply(tau_seq, function(th) {
    A <- ifelse(Tm >= th, 1, 0); diag(A) <- 0
    G <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
    if (gorder(G) == 0) return(tibble(tau = th, n_edges = 0, n_comp = 0, max_comp = 0, mean_deg = 0))
    comps <- components(G)
    tibble(tau = th,
           n_edges = gsize(G),
           n_comp  = length(comps$csize),
           max_comp = max(comps$csize),
           mean_deg = mean(degree(G)))
  }))
}

plot_threshold_sweep <- function(sweep_df, tag) {
  p1 <- ggplot(sweep_df, aes(tau, n_edges)) +
    geom_line() + geom_point() +
    labs(title = paste0("Edge count vs τ • ", tag), x = "τ threshold", y = "Edges") +
    theme_pub()
  p2 <- ggplot(sweep_df, aes(tau, max_comp)) +
    geom_line() + geom_point() +
    labs(title = paste0("Largest component size vs τ • ", tag), x = "τ threshold", y = "Largest component") +
    theme_pub()
  list(p1 = p1, p2 = p2)
}
# ======================================================


# ====== SUMMARY-STATS LINK (+ PS–PS τ/ρ/r HEATMAPS) ======
read_sumstats_and_cor <- function(path_to_csv, ps_cols, stat_start) {
  df <- read_csv(path_to_csv, show_col_types = FALSE)
  
  PS      <- df[, ps_cols,   drop = FALSE]
  sumstat <- df[, stat_start:ncol(df), drop = FALSE]
  
  # Method × Stat Kendall (for interpretability)
  ct_ms <- psych::corr.test(as.matrix(cbind(PS, sumstat)),
                            method = "spearman", adjust = "none")
  r_ms <- ct_ms$r[colnames(PS), colnames(sumstat)]
  p_ms <- ct_ms$p[colnames(PS), colnames(sumstat)]
  
  cor_df <- reshape2::melt(r_ms, varnames = c("PS","Stat"), value.name = "r") %>%
    left_join(
      reshape2::melt(p_ms, varnames = c("PS","Stat"), value.name = "p"),
      by = c("PS","Stat")
    ) %>%
    mutate(star = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ))
  
  # Method × Method: Kendall (τ), Spearman (ρ), Pearson (r)
  ct_tau <- psych::corr.test(as.matrix(PS), method = "kendall",  adjust = "none")
  ct_spr <- psych::corr.test(as.matrix(PS), method = "spearman", adjust = "none")
  ct_pea <- psych::corr.test(as.matrix(PS), method = "pearson",  adjust = "none")
  
  list(
    cor_df = cor_df,         # method × stat Kendall
    r_ms   = r_ms,
    # PS–PS correlations:
    r_mm_tau   = ct_tau$r,   p_mm_tau   = ct_tau$p,
    r_mm_spear = ct_spr$r,   p_mm_spear = ct_spr$p,
    r_mm_pear  = ct_pea$r,   p_mm_pear  = ct_pea$p,
    ps_names   = colnames(PS)
  )
}

plot_method_stat_heatmap_stars <- function(cor_df, method_order = NULL, tag = "") {
  df <- cor_df
  if (!is.null(method_order)) df <- df %>% mutate(PS = factor(PS, levels = method_order))
  ggplot(df, aes(x = Stat, y = PS, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_div_pub(lim = c(-1,1)) +
    geom_text(aes(label = sprintf("%.2f%s", r, star)), size = 3) +
    labs(
      title = paste0("Method × Summary-stat Kendall correlation (with p-stars) • ", tag),
      x = "Summary Statistic", y = "Plasticity Score / Method", fill = "τ"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
}

# Generic method×method correlation heatmap (used for τ/ρ/r)
plot_method_corr_heatmap <- function(M, order = NULL, tag = "", title = "", fill_label = "corr") {
  X <- M
  if (!is.null(order)) {
    keep <- intersect(order, rownames(X))
    X <- X[keep, keep, drop = FALSE]
  }
  df <- reshape2::melt(X, varnames = c("method_i","method_j"), value.name = "val")
  if (!is.null(order)) {
    ord <- order[order %in% unique(df$method_i)]
    df <- df %>%
      mutate(method_i = factor(method_i, levels = ord),
             method_j = factor(method_j, levels = ord))
  }
  ggplot(df, aes(method_j, method_i, fill = val)) +
    geom_tile(color = "white") +
    coord_fixed() +
    scale_fill_div_pub(lim = c(-1,1)) +
    labs(title = paste0(title, " • ", tag), x = NULL, y = NULL, fill = fill_label) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

# Build a method×method similarity from stat profiles (Pearson corr of r-vectors)
stat_profile_similarity_matrix_from_r <- function(r_ms) {
  if (length(r_ms) == 0) return(matrix(NA_real_, 0, 0))
  S <- suppressWarnings(cor(t(r_ms), use = "pairwise.complete.obs", method = "pearson"))
  dimnames(S) <- list(rownames(r_ms), rownames(r_ms))
  S
}

plot_matrix_link_scatter <- function(Tm, S, tag, x_lab = "Kendall’s τ", y_lab = "Stat-profile similarity") {
  common <- intersect(rownames(Tm), rownames(S))
  if (length(common) < 3) return(ggplot() + theme_void())
  
  T2 <- Tm[common, common, drop = FALSE]
  S2 <- S[common, common, drop = FALSE]
  lt <- lower.tri(T2)
  df <- tibble(tau = as.vector(T2[lt]), sim = as.vector(S2[lt])) %>%
    filter(!is.na(tau), !is.na(sim))
  
  r_spearman <- suppressWarnings(cor(df$tau, df$sim, method = "spearman", use = "complete.obs"))
  
  ggplot(df, aes(tau, sim)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE) +
    annotate("text",
             x = min(df$tau, na.rm = TRUE), y = max(df$sim, na.rm = TRUE),
             hjust = 0, vjust = 1,
             label = paste0("matrix ρ = ", sprintf("%.2f", r_spearman), " (Spearman)")) +
    labs(title = paste0("Linking τ to stat-profile similarity • ", tag),
         x = x_lab, y = y_lab) +
    theme_pub()
}

plot_method_stat_pca_biplot <- function(r_ms, Tm, tau_threshold, tag, alpha = 0.5) {
  # r_ms: methods × stats matrix (e.g., Kendall correlations of methods vs stats)
  if (nrow(r_ms) == 0) return(ggplot() + theme_void())
  
  # Center/scale stats across methods (PCA on correlation-like data)
  X <- scale(r_ms, center = TRUE, scale = TRUE)
  
  # SVD: X = U D V^T
  sv  <- svd(X)
  U   <- sv$u
  D   <- sv$d
  V   <- sv$v
  
  # Biplot scaling: F = U D^alpha, G = V D^(1-alpha)
  D_alpha   <- diag(D^alpha, nrow = length(D))
  D_1alpha  <- diag(D^(1 - alpha), nrow = length(D))
  F         <- U %*% D_alpha         # methods (rows of X)
  G         <- V %*% D_1alpha        # stats   (columns of X)
  
  # Grab first two PCs
  scores <- tibble(PC1 = F[,1], PC2 = F[,2], method = rownames(r_ms))
  loads  <- tibble(PC1 = G[,1], PC2 = G[,2], Stat   = colnames(r_ms))
  
  # Colour methods by τ-graph components (like before)
  A <- ifelse(Tm >= tau_threshold, 1, 0); diag(A) <- 0
  Gtau <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  comps <- if (igraph::gorder(Gtau) > 0) igraph::components(Gtau)$membership else
    setNames(rep(1, nrow(r_ms)), rownames(r_ms))
  scores <- scores %>% mutate(component = factor(unname(comps[method])))
  
  # Put arrows on a comparable scale to points
  s_rng <- max(abs(dplyr::select(scores, PC1, PC2)))
  l_rng <- max(abs(dplyr::select(loads,  PC1, PC2)))
  arrow_scale <- if (l_rng > 0) 0.9 * s_rng / l_rng else 1
  loads <- loads %>% mutate(PC1 = PC1 * arrow_scale, PC2 = PC2 * arrow_scale)
  
  ggplot() +
    # statistic arrows
    geom_segment(data = loads,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.12, "inches")), alpha = 0.6) +
    ggrepel::geom_text_repel(data = loads, aes(x = PC1, y = PC2, label = Stat),
                             size = 3, max.overlaps = 50) +
    # method points
    geom_point(data = scores, aes(x = PC1, y = PC2, fill = component),
               shape = 21, size = 3, color = "black") +
    ggrepel::geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = method, color = component),
                             size = 3, max.overlaps = 50, show.legend = FALSE) +
    scale_fill_manual(values = scales::hue_pal(l = 55)(10), name = "τ-component") +
    scale_color_manual(values = scales::hue_pal(l = 55)(10), guide = "none") +
    labs(
      title = paste0("PCA biplot of method–stat relationships • ", tag,
                     " (α = ", alpha, ")"),
      x = "PC1", y = "PC2"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# =======================================================================


# ===================== PIPELINE LOOP ===================
pairs_files <- find_files(in_dir, "^kendall_pairs_.*_res\\d+\\.csv$")
rank_files  <- find_files(in_dir, "^ranks_matrix_combined_.*_res\\d+\\.csv$")
if (length(pairs_files) == 0 || length(rank_files) == 0) stop("Missing inputs in ", in_dir)

pdf_path <- file.path(out_dir, "analysis_report.pdf")
pdf(pdf_path, width = 10, height = 7, onefile = TRUE)

for (pairs_path in pairs_files) {
  tag <- build_tag(pairs_path, "kendall_pairs")
  ranks_path <- rank_files[grepl(tag, rank_files)]
  if (length(ranks_path) == 0) next
  ranks_path <- ranks_path[1]
  
  grid::grid.newpage()
  grid::grid.text(paste0("Analysis report • ", tag), x = 0.5, y = 0.9,
                  gp = grid::gpar(cex = 1.6, fontface = "bold"))
  grid::grid.text(paste0("τ threshold = ", tau_threshold,
                         "    k = {", paste(topk_values, collapse = ", "), "}"),
                  x = 0.5, y = 0.82, gp = grid::gpar(cex = 1.0))
  
  # --- Kendall pairs -> matrix
  pairs_df <- read_csv(pairs_path, show_col_types = FALSE) |>
    mutate(correlation = as.numeric(correlation)) |>
    filter(!is.na(correlation))
  idx <- sort(unique(c(pairs_df$score_i, pairs_df$score_j)))
  Tm <- matrix(NA_real_, length(idx), length(idx), dimnames = list(idx, idx))
  diag(Tm) <- 1
  for (k in seq_len(nrow(pairs_df))) {
    i <- pairs_df$score_i[k]; j <- pairs_df$score_j[k]; v <- pairs_df$correlation[k]
    Tm[i, j] <- v; Tm[j, i] <- v
  }
  
  # --- τ visuals
  print(plot_tau_heatmap(Tm, tag))
  sw <- threshold_sweep_stats(Tm, tau_seq = seq(0.3, 0.9, by = 0.05))
  sw_plots <- plot_threshold_sweep(sw, tag); print(sw_plots$p1); print(sw_plots$p2)
  print(plot_agreement_network(Tm, tau_threshold, tag))
  print(plot_component_sizes(Tm, tau_threshold, tag))
  
  # --- Ranks
  ranks_df <- read_csv(ranks_path, show_col_types = FALSE) |> ensure_gid()
  rank_cols <- setdiff(colnames(ranks_df), c("gid","genotype_id"))
  ranks_df <- ranks_df |> mutate(across(all_of(rank_cols), as.numeric))
  
  # --- Flips and top-k
  flips_tbl <- if (length(rank_cols) >= 2) {
    combn(rank_cols, 2, simplify = FALSE) |>
      purrr::map_dfr(function(pr) {
        a <- pr[1]; b <- pr[2]
        sub_df <- dplyr::select(ranks_df, dplyr::all_of(c(a, b)))
        sub_df <- sub_df[stats::complete.cases(sub_df), , drop = FALSE]
        ef <- exact_flips_from_ranks(sub_df[[a]], sub_df[[b]])
        tibble(tag = tag, score_i = a, score_j = b,
               flips_abs = ef$discordant, n_pairs = ef$n_pairs, flips_frac = ef$flips_frac)
      })
  } else tibble(tag = character(), score_i = character(), score_j = character(),
                flips_abs = integer(), n_pairs = integer(), flips_frac = numeric())
  write_csv(flips_tbl, file.path(out_dir, paste0("rank_flips_", tag, ".csv")))
  if (nrow(flips_tbl) > 0) print(plot_flips_vs_tau(flips_tbl, pairs_df, tag))
  
  topk_df <- topk_overlap_tbl(ranks_df, topk_values) |> mutate(tag = tag, .before = 1)
  write_csv(topk_df, file.path(out_dir, paste0("topk_overlap_", tag, ".csv")))
  
  metrics <- method_agreement_metrics(Tm, tau_threshold, topk_df, flips_tbl)
  write_csv(metrics, file.path(out_dir, paste0("method_agreement_metrics_", tag, ".csv")))
  print(plot_method_agreement(metrics, tau_threshold, tag))
  print(plot_method_flips(metrics, tag))
  
  mos <- method_overlap_summary_from_topk(topk_df, agree_thr = AGREE_THR_FULL)
  write_csv(mos, file.path(out_dir, paste0("method_overlap_summary_global_", tag, ".csv")))
  print(plot_method_k_heatmap(mos, tag, value_col = "avg_overlap",
                              title_suffix = "Average top-k overlap vs others (global)"))
  print(plot_method_k_heatmap(mos, tag, value_col = "min_overlap",
                              title_suffix = "Minimum top-k overlap vs others (global)"))
  
  k_show <- max(topk_values)
  Om <- topk_overlap_matrix(ranks_df, rank_cols, k_show)
  print(plot_topk_overlap_heatmap_matrix(Om, tag, k_show))
  
  # ================= PS–STAT LINK + PS–PS CORR HEATMAPS =================
  if (!is.null(sumstats_path_override) && file.exists(sumstats_path_override)) {
    ss <- read_sumstats_and_cor(sumstats_path_override,
                                ps_cols = sumstats_ps_cols,
                                stat_start = sumstats_stat_from)
    cor_df    <- ss$cor_df
    r_ms      <- ss$r_ms
    r_mm_tau  <- ss$r_mm_tau
    r_mm_spr  <- ss$r_mm_spear
    r_mm_pear <- ss$r_mm_pear
    
    # Use τ heatmap order for comparability
    m <- Tm; diag(m) <- 1; m[is.na(m)] <- 0
    ord <- { d <- as.dist(1 - m); hc <- hclust(d, method = "average"); rownames(m)[hc$order] }
    
    # Keep only methods common across sources
    common_methods <- Reduce(intersect, list(ord, rownames(r_ms), rownames(r_mm_tau)))
    cor_df    <- cor_df %>% filter(PS %in% common_methods)
    r_ms      <- r_ms[common_methods, , drop = FALSE]
    r_mm_tau  <- r_mm_tau[common_methods, common_methods, drop = FALSE]
    r_mm_spr  <- r_mm_spr[common_methods, common_methods, drop = FALSE]
    r_mm_pear <- r_mm_pear[common_methods, common_methods, drop = FALSE]
    ord       <- ord[ord %in% common_methods]
    
    # (1) Method × stat Kendall with stars
    print(plot_method_stat_heatmap_stars(cor_df, method_order = ord, tag = tag))
    # (2) Method × method Kendall (PS realisations) — as before
    print(plot_method_corr_heatmap(r_mm_tau, order = ord, tag = tag,
                                   title = "Kendall correlation (τ) between plasticity score realisations",
                                   fill_label = "τ"))
    
    # (3) Link τ to stat-profile similarity
    S <- stat_profile_similarity_matrix_from_r(r_ms)
    write_csv(as_tibble(S, rownames = "method"),
              file.path(out_dir, paste0("stat_profile_similarity_matrix_", tag, ".csv")))
    print(plot_matrix_link_scatter(Tm[common_methods, common_methods], S, tag))
    
    # (4) PCA biplot of method–stat Kendall correlations
    print(plot_method_stat_pca_biplot(r_ms, Tm[common_methods, common_methods], tau_threshold, tag))
    
    # === NEW PAGES AT THE END: Spearman & Pearson PS–PS HEATMAPS ===
    print(plot_method_corr_heatmap(r_mm_spr, order = ord, tag = tag,
                                   title = "Spearman correlation (ρ) between plasticity score realisations",
                                   fill_label = "ρ"))
    print(plot_method_corr_heatmap(r_mm_pear, order = ord, tag = tag,
                                   title = "Pearson correlation (r) between plasticity score realisations",
                                   fill_label = "r"))
  } else {
    grid::grid.newpage()
    grid::grid.text("Summary-statistics file not found.\nSet 'sumstats_path_override' to your CSV.",
                    x = 0.5, y = 0.5, gp = grid::gpar(cex = 1))
  }
  # ==================================================================
  
  cat("Completed:", tag, "\n\n")
}



dev.off()
message("Wrote PDF: ", pdf_path)
writeLines("Done.")
