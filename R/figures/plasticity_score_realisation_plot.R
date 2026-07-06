# =============================================================================
# plasticity_score_realisation_plot.R — plasticity-score realisations (simulated)
# =============================================================================
# WHAT IT DOES: Plots the realised plasticity scores across genotypes/forms for the
#   simulated data.
# REQUIRES:     Score outputs produced by 03_plasticity_scores.R (run a scenario first),
#               read from under output/.
# PRODUCES:     Plasticity-score realisation figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plasticity_score_realisation_plot.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named input path / selection assignments in the body)
#   the input score path and any index/form selection are set in the script body.
# =============================================================================


#this script is for the assembly of the second figure containing the boxplots for the distribution of the scores as well as the clustering results (hclust with kendalll correlation distance as distance metrics) in form of a dendrogram 
#and the adjusted rand index table comparing the clusterings across the different genotype forms. THe reszlts suggest (sampling interval=5) that the scores obey to a similar clustering for all different genotype forms

library(dplyr)
library(readr)
library(cluster)      # silhouette()
library(mclust)       # adjustedRandIndex()
library(factoextra)   # fviz_dend()
library(patchwork)    # layou
library(cluster)     # for silhouette()
library(pheatmap)
library(dendextend)  # for nicer dendrogram handling
library(ggbeeswarm)
library(tidyr)
library(ggplot2)
library(patchwork)
library(forcats)    # for fct_reorder()
library(ggridges)   # for geom_density_ridges()
library(pheatmap)   # for heatmaps
library(ggpubr)    # for ggtexttable & ggarrange
library(dplyr)     # for mutate_if
library(tibble)    # for rownames_to_column
library(aricode)
library(knitr)
library(scales)
library(vegan)      # for mantel()
# Read 

if (exists("form_ranges") && length(form_ranges) == 3) {
  scores_list <- readRDS(here::here("data", "maize_scores_list.rds"))
  scores_list_arr = array(unlist(scores_list),dim=c(732,3,28))
  scores_list_arr2 = as.numeric(scores_list_arr[,1,])
  dim(scores_list_arr2) = c(732,28)
  colnames(scores_list_arr2) = names(scores_list)
  df = data.frame(scores_list_arr2)
  df$Genotype = seq(732)
} else {
  form_ranges <- list(
    linear     = 1:20,
    gaussian   = 21:40,
    sinusoidal = 41:60,
    wave       = 61:80
  )

  df <- read_csv(here::here("data", "regression_summary_stats", "regression_data_full_interval_10_indices_1_11_21_31_41_50.csv"))
  na_counts <- colSums(is.na(df))

  # keep only columns with zero NAs (there are some in the case of 1 and 2 samples across the environmental gradient as the RN and RNN cannot be calculated in this scenario)
  df <- df[ , na_counts == 0]

  df <- df %>% rename(Genotype = names(df)[1])
  stat_cols <- c("min","max","mean","median","slope","range","variance","mean_lower","mean_upper")
  df <- df[, !colnames(df) %in% stat_cols, drop = FALSE]
}

scores_long <- df %>%
  pivot_longer(
    cols       = -Genotype,       # all the metric columns
    names_to   = "Metric",
    values_to  = "Score"
  )

# Make one histogram + density per metric. This is a huge figure alone with all the distributions of all the plasticity scores. I suppose this can either go into the supplements or be left out entirely
make_plot <- function(metric_name) {
  sub <- scores_long %>% filter(Metric == metric_name)
  ggplot(sub, aes(x = Score)) +
    geom_histogram(aes(y = ..density..),
                   bins  = 30,
                   fill  = "grey70",
                   color = "white") +
    geom_density(size = 1, color = "black") +
    labs(x = NULL, y = "Density") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
}


plots   <- lapply(unique(scores_long$Metric), make_plot)
combined <- wrap_plots(plots, ncol = 4) 
print(combined)
ggsave(create.dir = TRUE, here::here("output", "plots", "plasticity_scores_histograms.pdf"), combined, width = 6.3, height = 5,
       dpi = 900, units = "in", device = "pdf")


##############################
scores_long_trimmed <- scores_long %>%
  group_by(Metric) %>%
  filter(Score <= quantile(Score, 0.95, na.rm = TRUE)) %>%
  ungroup()

pp=ggplot(scores_long_trimmed, aes(x = Metric, y = Score)) +
  geom_boxplot(
    width        = 0.9,
    outlier.size = 0.5,
    fill         = "grey80",
    color        = "grey30"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = NULL,
    y = "Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.length  = unit(1, "pt"),
    plot.margin        = margin(2, 2, 2, 2)
  )
print(pp)

ggsave(create.dir = TRUE, here::here("output", "plots", "figures", "plasticity_scores_boxplots.pdf"), pp, width = 6.3, height = 5,
       dpi = 900, units = "in", device = "pdf")

##########################################################


# Zoran suggested that this might look good with violin plots but i think it doesnt 
pp_violin = ggplot(scores_long_trimmed, aes(x = Metric, y = Score)) +
  geom_violin(
    trim   = FALSE,
    width  = 1.2,
    adjust = 2,
    fill   = "grey80",
    color  = "grey30"
  ) +
  geom_jitter(
    aes(color = Metric),
    position    = position_jitter(width = 0.2, height = 0),
    size        = 0.7,
    alpha       = 0.5,
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = NULL,
    y = "Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.length  = unit(1, "pt"),
    plot.margin        = margin(2, 2, 2, 2)
  )

print(pp_violin)
############################################################

df=df[,!names(df)=="Genotype"]





forms      <- names(form_ranges)
all_labels <- c(forms, "combined")

# --- Build hclust objects and distance matrices for each form ---
hc_list   <- list()
dist_list <- list()

for (f in forms) {
  sub_df       <- df[ form_ranges[[f]], , drop = FALSE ]
  C            <- cor(sub_df, method = "kendall", use = "pairwise.complete.obs")
  D            <- as.dist(1 - C)
  dist_list[[f]] <- D
  hc_list[[f]]   <- hclust(D, method = "complete")
}

C_all                  <- cor(df, method = "kendall", use = "pairwise.complete.obs")
D_all                  <- as.dist(1 - C_all)
dist_list[["combined"]] <- D_all
hc_list[["combined"]]   <- hclust(D_all, method = "complete")

# --- Select k via average silhouette width across all four forms ---
k_range <- 3:10
sil_profile <- sapply(k_range, function(ki) {
  mean(sapply(forms, function(f) {
    cl  <- cutree(hc_list[[f]], k = ki)
    mean(silhouette(cl, dist_list[[f]])[, "sil_width"])
  }))
})

k <- k_range[which.max(sil_profile)]
message("Optimal k by average silhouette width: ", k)

# --- Silhouette profile plot (diagnostic) ---
sil_df <- data.frame(k = k_range, avg_sil = sil_profile)
p_sil <- ggplot(sil_df, aes(x = k, y = avg_sil)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = k_range) +
  labs(x = "Number of clusters (k)", y = "Average silhouette width") +
  theme_minimal(base_size = 10)
print(p_sil)

# --- Cut all dendrograms at optimal k (for dendrogram colouring only) ---
clusterings <- list()
for (lb in all_labels) {
  clusterings[[lb]] <- cutree(hc_list[[lb]], k = k)
}

# --- Mantel test: compare 27x27 score distance matrices across forms ---
# Tests whether pairwise inter-score distances are correlated between forms,
# without requiring any choice of k.
nA       <- length(all_labels)
mantel_r <- matrix(NA, nA, nA, dimnames = list(all_labels, all_labels))
mantel_p <- mantel_r

for (i in seq_len(nA - 1)) {
  for (j in (i + 1):nA) {
    a   <- all_labels[i]
    b   <- all_labels[j]
    res <- mantel(dist_list[[a]], dist_list[[b]], method = "kendall", permutations = 999)
    mantel_r[a, b] <- mantel_r[b, a] <- res$statistic
    mantel_p[a, b] <- mantel_p[b, a] <- res$signif
  }
}


# --- Manual functional grouping of the scores (NOT the hclust clusters) ---
# Shown as a shape marker to the RIGHT of each (cluster-coloured) PPI label:
#   Dispersion -> cross, Slope -> square, Range -> circle.
score_groups <- list(
  "Dispersion" = c("CV_t","CEV","SI","RSI","gPi","ESP","PFI","RDPI","RPI","EVS","PPF"),
  "Slope"      = c("RN","FW","PSI","RNN","ESPIID","RTR","APC","PQ","ESPI"),
  "Range"      = c("PR","NRW","PPi","PImd","PILSM","PIR","D_slope","RC")
)

score_group <- setNames(
  rep(names(score_groups), lengths(score_groups)),
  unlist(score_groups)
)
group_levels <- names(score_groups)
group_shapes <- c(Dispersion = 4, Slope = 0, Range = 16)   # cross / hollow square / circle

make_dend_plot <- function(hc, title, shrink = 0.5) {
  # Shrinks AND reverses the axis to keep the root on the left
  shrink_rev_trans <- trans_new(
    name      = paste0("shrink_rev_", shrink),
    transform = function(x) -(x * shrink),
    inverse   = function(x) -(x / shrink)
  )
  
  p <- fviz_dend(hc,
                 k           = k,
                 rect        = FALSE,
                 show_labels = TRUE,   # cluster-coloured PPI labels (k_colors)
                 repel       = FALSE,
                 k_colors    = "jco",
                 cex         = 0.6,
                 main        = title,
                 horiz       = TRUE    # <-- CRITICAL: Forces the text labels to be horizontal (angle = 0)
  )

  # Functional-group shape marker, placed in a column just to the right of the
  # labels. The shrink_rev transform runs Height left->right as 2 .. 0 .. -1,
  # so the label gutter is at negative y; marker_y clears the longest label.
  marker_y <- -0.95
  leaf_df <- data.frame(
    x   = seq_along(hc$order),
    lab = hc$labels[hc$order],
    stringsAsFactors = FALSE
  )
  leaf_df$grp <- factor(score_group[leaf_df$lab], levels = group_levels)

  p <- p +
    geom_point(
      data        = leaf_df,
      mapping     = aes(x = x, y = marker_y, shape = grp),
      inherit.aes = FALSE,
      colour      = "black",
      size        = 1.6
    ) +
    scale_shape_manual(values = group_shapes, name = "Functional group",
                       drop = FALSE,
                       guide = guide_legend(nrow = 1)) +

    scale_y_continuous(trans = shrink_rev_trans) +

    coord_flip(clip = "off") +
    
    theme_minimal(base_size = 10) +
    theme(
      aspect.ratio = 1.5,
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      
      # --- CRITICAL AXIS FIXES ---
      axis.ticks.y = element_blank(),  # Hide leaf index ticks on the left
      axis.text.y  = element_blank(),  # Hide leaf index labels (0, 10, 20) on the left
      
      axis.ticks.x = element_blank(),   # Show height ticks on the bottom (optional but recommended)
      axis.text.x  = element_blank(),   # Show height labels on the bottom (visual x-axis)
      # ---------------------------
      
      # Right margin accommodates the long horizontal labels + shape markers
      plot.margin  = margin(t = 10, r = 28, b = 5, l = 0),

      # Group-membership legend at the bottom (collected once across panels)
      legend.position   = "bottom",
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7),
      legend.key.size   = unit(3, "mm")
    )

  return(p)
}



plots_sep <- lapply(forms, function(f) {
  make_dend_plot(hc_list[[f]], title = f)
})
plot_comb <- make_dend_plot(hc_list[["combined"]], title = "combined")
print(plot_comb)

if (length(plots_sep) == 3) {
sep_panel <- (plot_comb | plots_sep[[2]]) /
  (plots_sep[[3]])
} else {
sep_panel <- (plot_comb | plots_sep[[2]]) /
  (plots_sep[[3]] | plots_sep[[4]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

print(sep_panel)
ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "plasticity_scores_dendrograms.pdf"),
  sep_panel,
  width = 6.3, height = 5, dpi = 900, units = "in", device = "pdf")

ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "dendrogram_linear.pdf"),
  plots_sep[[1]],
  width = 6.3, height = 4, dpi = 900, units = "in", device = "pdf")

if (length(plots_sep) == 4) {
form_panel <- (plots_sep[[1]] | plots_sep[[2]]) /
  (plots_sep[[3]] | plots_sep[[4]])
}

########################################


sig_asterisks <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

obs_df <- mantel_r %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form")

p_df <- mantel_p %>%
  as.data.frame() %>%
  rownames_to_column(var = "Form")

obs_annot_df <- obs_df
for (col in names(obs_annot_df)[-1]) {
  obs_annot_df[[col]] <- sprintf(
    "%.3f%s",
    obs_df[[col]],
    sig_asterisks(p_df[[col]])
  )
}

obs_long <- obs_annot_df %>%
  pivot_longer(-Form, names_to = "Comparison", values_to = "Mantel_r_star")


ord <- all_labels 

tbl_obs_plot <- ggplot(obs_long, aes(x = Comparison, y = Form)) +
  scale_x_discrete(limits = ord, position = "top") +
  scale_y_discrete(limits = ord) +

  geom_tile(aes(fill = as.numeric(sub("\\*+$", "", Mantel_r_star))), color = "grey70") +
  scale_fill_distiller(palette = "OrRd", direction = 1,
                       limits = c(0, 1), na.value = "grey85",
                       name = "Mantel r") +


  geom_text(aes(label = gsub("[^\\*]", "", Mantel_r_star)), size = 5, vjust = 0.75) +


coord_fixed() +

  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    direction      = "horizontal",
    barwidth       = unit(40, "mm"),
    barheight      = unit(3, "mm")
  )) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),

    axis.text.x     = element_text(angle = 45, hjust = 0, vjust = 0, face = "bold"),
    axis.text.y     = element_text(face = "bold"),
    axis.title      = element_blank(),
    legend.position = "bottom",
    plot.margin     = margin(2, 2, 2, 2, "mm")
  )
print(tbl_obs_plot)

ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "plasticity_scores_mantel_table.pdf"),
  tbl_obs_plot,
  width = 6.3, height = 5, dpi = 900, units = "in", device = "pdf"
)


####################




final_fig3 <- wrap_plots(
  A = wrap_elements(full = pp),
  B = wrap_elements(full = plots_sep[[1]]),
  C = wrap_elements(full = tbl_obs_plot),
  design = "
AAAAAA
AAAAAA
AAAAAA
BBBCCC
BBBCCC
BBBCCC
"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )


print(final_fig3)

ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "figure2_2.pdf"),
  final_fig3,
  width = 6.3, height = 8, dpi = 900, units = "in", device = "pdf"
)


################



final_fig4 <- wrap_plots(
  A = wrap_elements(full = pp_violin),
  B = wrap_elements(full = plots_sep[[1]]),
  #C = wrap_elements(full = plot_comb),
  C= wrap_elements(full = tbl_obs_plot),
  design = "
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBBBCC
BBBBCC

"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "figure2_3.pdf"),
  final_fig4,
  width = 6.3, height = 8, dpi = 900, units = "in", device = "pdf"
)
print(final_fig4)
#########################################


#this figure is supposed to go in the supplements
if (length(plots_sep) == 4) {
final_fig7 <- wrap_plots(
  A = wrap_elements(full = plots_sep[[2]]),
  B = wrap_elements(full = plots_sep[[3]]),
  C= wrap_elements(full = plots_sep[[4]]),
  D = wrap_elements(full = plot_comb),
  design = "
AABB
AABB
CCDD
CCDD
"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )
} else {
final_fig7 <- wrap_plots(
  A = wrap_elements(full = plots_sep[[2]]),
  B = wrap_elements(full = plots_sep[[3]]),
  C = wrap_elements(full = plot_comb),
  design = "
AABB
AABB
CCDD
CCDD
"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )
}

ggsave(create.dir = TRUE, 
  here::here("output", "plots", "figures", "dendro_supp.pdf"),
  final_fig7,
  width = 6.3, height = 6, dpi = 900, units = "in", device = "pdf"
)

print(final_fig7)
