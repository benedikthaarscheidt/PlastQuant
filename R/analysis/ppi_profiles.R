# =============================================================================
# ppi_profiles.R — groups of PPIs by what they measure (correlation profiles)
# =============================================================================
# Every PPI is described by its profile: its signed correlations with properties of the
# reaction norms (summary statistics, or the true parameters of a simulation). PPIs with
# similar profiles measure the same properties of the reaction norm. These groups are
# compared with the groups of PPIs that rank the genotypes alike (Kendall clustering, as
# in the dendrogram of plasticity_score_realisation_plot.R).

# Hierarchical clustering of the PPIs on their profiles. P: matrix, rows = PPIs, columns =
# properties, values = signed correlations; properties without a value for every PPI are
# left out. Euclidean distance, average linkage; number of groups in 2 ... max_groups with
# the largest average silhouette width. Returns the distance, tree, number and groups.
profile_clustering <- function(P, max_groups = 8) {
  P <- P[, colSums(is.finite(P)) == nrow(P), drop = FALSE]
  d <- dist(P)
  hc <- hclust(d, method = "average")
  ks <- 2:min(max_groups, nrow(P) - 1)
  sil <- vapply(ks, function(g) mean(cluster::silhouette(cutree(hc, g), d)[, "sil_width"]), numeric(1))
  k <- ks[which.max(sil)]
  list(profile = P, dist = d, hc = hc, k = k, groups = cutree(hc, k))
}

# Profile of every PPI from summary statistics: Spearman correlation over all rows of tab
# (one row per genotype, e.g. population_score_table()) between every PPI and every
# statistic. Rows of the result: ppis, columns: stats.
summary_stat_profile <- function(tab, ppis, stats = c("mean", "min", "max", "range", "slope")) {
  sapply(stats, function(s) vapply(ppis, function(p) suppressWarnings(cor(tab[[p]], tab[[s]],
         method = "spearman", use = "pairwise.complete.obs")), numeric(1)))
}

# Groups of PPIs that rank the genotypes alike: Kendall distance 1 - tau over all genotypes
# pooled, complete linkage, number of groups in `groups` with the largest average silhouette
# width averaged over the per-shape clusterings. df: genotypes x PPIs; rows_by_shape: list
# of the rows of every shape.
kendall_ppi_clustering <- function(df, rows_by_shape, groups = 3:10) {
  kd <- function(x) as.dist(1 - cor(x, method = "kendall", use = "pairwise.complete.obs"))
  D <- lapply(rows_by_shape, function(r) kd(df[r, , drop = FALSE]))
  hcs <- lapply(D, hclust, method = "complete")
  sil <- vapply(groups, function(g) mean(vapply(names(D), function(s)
    mean(cluster::silhouette(cutree(hcs[[s]], g), D[[s]])[, "sil_width"]), numeric(1))), numeric(1))
  C <- groups[which.max(sil)]
  d_all <- kd(df)
  hc <- hclust(d_all, method = "complete")
  list(dist = d_all, hc = hc, k = C, groups = cutree(hc, C))
}

# Colours of the groups as in the dendrograms (factoextra::fviz_dend, palette "jco"): the
# groups are coloured in the order in which they appear along the tree. Named by PPI.
group_colours <- function(hc, groups) {
  ord <- unique(groups[hc$order])
  pal <- setNames(ggpubr::get_palette("jco", length(ord)), ord)
  setNames(unname(pal[as.character(groups)]), names(groups))
}

# Agreement of the profile groups with the Kendall groups: Mantel test between the two
# distance matrices (Kendall correlation, 999 permutations) and adjusted Rand index of the
# two group assignments, over the PPIs present in both.
compare_ppi_clusterings <- function(prof, kend) {
  ppis <- intersect(labels(prof$dist), labels(kend$dist))
  sub <- function(d) as.dist(as.matrix(d)[ppis, ppis])
  m <- vegan::mantel(sub(prof$dist), sub(kend$dist), method = "kendall", permutations = 999)
  list(mantel_r = unname(m$statistic), mantel_p = m$signif,
       ari = aricode::ARI(prof$groups[ppis], kend$groups[ppis]))
}

# Figure: dendrogram of the profile clustering, root on the left, PPI names coloured by
# profile group (as the dendrograms of plasticity_score_realisation_plot.R); a square after
# every name shows the PPI's Kendall group (kend_col: colour per PPI, named by PPI).
# cmp: compare_ppi_clusterings(), reported in the subtitle.
plot_profile_clustering <- function(prof, kend_col, cmp, title, shrink = 0.5) {
  hc <- prof$hc
  H <- max(hc$height)
  shrink_rev <- scales::trans_new(paste0("shrink_rev_", shrink),
                                  transform = function(x) -(x * shrink),
                                  inverse = function(x) -(x / shrink))
  leaves <- hc$labels[hc$order]
  marks <- data.frame(x = seq_along(leaves), y = -0.3 * H, col = unname(kend_col[leaves]))
  factoextra::fviz_dend(hc, k = prof$k, k_colors = "jco", horiz = TRUE, cex = 0.6,
                        rect = FALSE, main = title) +
    ggplot2::geom_point(data = marks, ggplot2::aes(x = x, y = y, fill = col), inherit.aes = FALSE,
                        shape = 22, size = 2.2, colour = "grey30", stroke = 0.2) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(trans = shrink_rev) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::labs(subtitle = sprintf("%d profile groups (name colours) vs %d Kendall groups (squares);\nvs Kendall: Mantel r = %.2f, p = %.3f; adjusted Rand index = %.2f",
                                     prof$k, length(unique(kend_col)), cmp$mantel_r, cmp$mantel_p, cmp$ari)) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 10), plot.subtitle = ggplot2::element_text(size = 7),
                   plot.margin = ggplot2::margin(5, 25, 5, 5),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}
