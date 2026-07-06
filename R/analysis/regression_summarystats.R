# =============================================================================
# regression_summarystats.R — regress plasticity scores on summary statistics
# =============================================================================
# WHAT IT DOES: For each plasticity score, fits full/simple/ridge regression models on
#   the nine summary-stat predictors, builds coefficient tables with significance stars
#   and R^2, clusters scores, and writes per-score and combined diagnostic PDFs.
# REQUIRES:     regression_sumstats_dataprep.R MUST be run first — this script reads its
#               combined CSV from output/regression_summary_stats/.
# PRODUCES:     Coefficient tables and diagnostic PDFs under output/plots/.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","analysis","regression_summarystats.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   input_folder     path, default output/regression_summary_stats  where the CSV is read [COMMON]
#   output_folder    path, default output/plots                     where PDFs are written [COMMON]
#   predictor_names  character vector  the summary-stat predictors used
#   k                integer, default 4     number of score clusters
#   n_perm           integer, default 1000  permutations for the cluster test
# =============================================================================

# Script:    regression_summarystats.R
# Purpose:   Using a pre‑computed table of summary statistics versus plasticity scores,
#            this script fits and evaluates multiple regression models, generates
#            coefficient tables, diagnostic plots, and combined summary PDFs.
#
# For each score (column in the input CSV):
#   1. Split data 80/20 into training and test sets.
#   2. Fit a “full” linear model with all nine summary‑stat predictors, dropping any
#      collinear terms.
#   3. Fit three “simple”/specialized models:  
#        • univariate model with `range` only  
#        • univariate model with `mean_upper` only  
#        • ridge regression with small λ  
#   4. Extract and “tidy” model coefficients, append significance stars (*,**,***),
#      and compute R² & adjusted‑R².
#   5. Produce a per‑score PDF containing:  
#        • Table of full‑model estimates + R² & Adj‑R² rows  
#        • Table of simple‑model estimates + stars  
#        • Table of simple‑model adjusted R²  
#        • Predicted vs. observed scatter plots for all four models  
#        • Coefficient‑estimate plots with 95% CIs  
#        • Visreg effect plots for each predictor in the full model  
#        • Observed vs. fitted lines for the two best single‑predictor models  
#
# After looping all scores:
#   • Assemble “wide” tables of full‑model estimates and simple‑model estimates
#     (terms × scores), with stars.
#   • Build matrices of full‑model R² & Adj‑R² (scores × {R²,AdjR²}).
#   • Build matrix of simple‑model adjusted R² (predictors × scores).
#   • Generate a barplot of average simple‑model R² per predictor.
#   • Perform hierarchical clustering of scores based on their simple‑model R²
#     profiles, plot a dendrogram, compute silhouette widths, and run permutation
#     tests on cluster and overall silhouette means.
#   • Save all combined tables, plots, and clustering output into a multi‑page PDF.
#
# Inputs:
#   • CSV (e.g. regression_data_full_interval_2_…csv) with one row per observation,
#     columns = summary stats + plasticity scores.
#   • `input_folder` / `output_folder` paths.
# Dependencies:
#   MASS, car, ggplot2, broom, gridExtra, grid, visreg, dplyr, purrr, tidyr,
#   mclust, cluster, dendextend, RColorBrewer
# Outputs:
#   • One PDF per score (full/simple tables + plots) in `output_folder`
#   • One combined PDF (“all_scores_wide.pdf”) with coefficient tables, R² tables,
#     average R² barplot, dendrogram, silhouette tables, and permutation results.
# Load required libraries
library(MASS)
library(car)
library(ggplot2)
library(broom)
library(gridExtra)
library(grid)
library(visreg)
library(dplyr)
library(purrr)
library(tidyr)
library(mclust)
library(cluster)
library(dendextend)
library(RColorBrewer)

ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

# Directories for input data and output
input_folder  <- here::here("data", "regression_summary_stats")
output_folder <- here::here("output", "plots")

ensure_dir_exists(input_folder)
ensure_dir_exists(output_folder)

# Load combined data
combined_file <- file.path(input_folder,
                           "regression_data_full_interval_2_indices_1_3_5_7_9_11_13_15_17_19_21_23_25_27_29_31_33_35_37_39_41_43_45_47_49_50.csv")
combined_df <- read.csv(combined_file, row.names = 1)

# Predictor and score names
predictor_names <- c("min","max","mean","median","slope","range",
                     "variance","mean_lower","mean_upper")
score_names     <- setdiff(names(combined_df), predictor_names)

# Vectorized significance‐stars helper
signif_stars <- function(p) {
  sapply(p, function(x) {
    if      (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else                ""
  }, USE.NAMES = FALSE)
}

# Helper: Predicted vs Observed plot
plot_pred_obs <- function(obs, pred, title) {
  ggplot(data.frame(obs, pred), aes(x = pred, y = obs)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Predicted", y = "Observed", title = title) +
    theme_minimal()
}

# Helper: Coefficient plot
plot_coefs <- function(df, title) {
  df2 <- subset(df, term != "(Intercept)")
  df2$term <- factor(df2$term, levels = df2$term[order(df2$estimate)])
  p <- ggplot(df2, aes(x = term, y = estimate)) +
    geom_point(size = 3) + coord_flip() +
    labs(x = NULL, y = "Estimate", title = title) +
    theme_minimal()
  if ("std.error" %in% names(df2)) {
    p <- p + geom_errorbar(aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ), width = 0.2)
  }
  p
}

# Theme for tables
tt_small <- ttheme_minimal(
  core    = list(fg_params = list(fontsize = 6)),
  colhead = list(fg_params = list(fontsize = 7, fontface = "bold"))
)

set.seed(123)
full_models <- list()

# Per-score processing
for (score in score_names) {
  message("Processing score: ", score)
  
  # Data split
  df <- data.frame(y = combined_df[[score]], combined_df[, predictor_names])
  train_i <- sample(seq_len(nrow(df)), floor(0.8 * nrow(df)))
  df_tr   <- df[train_i, ]
  df_te   <- df[-train_i, ]
  
  # Full model with collinearity handling
  tmp_mod <- lm(y ~ ., data = df_tr)
  keep    <- setdiff(predictor_names, names(which(is.na(coef(tmp_mod)))))
  full_mod<- lm(as.formula(paste("y ~", paste(keep, collapse = " + "))), data = df_tr)
  
  # Other models
  mod_range<- lm(y ~ range, data = df_tr)
  mod_mu   <- lm(y ~ mean_upper, data = df_tr)
  mod_ridge<- lm.ridge(y ~ ., data = df_tr, lambda = 0.001)
  
  full_models[[score]] <- full_mod
  
  # Tidy outputs and compute R²
  td_full  <- tidy(full_mod)
  td_range <- tidy(mod_range)
  td_mu    <- tidy(mod_mu)
  td_ridge <- tidy(mod_ridge)
  s <- summary(full_mod)
  r2    <- round(s$r.squared,    2)
  adjr2 <- round(s$adj.r.squared, 2)
  
  # Page 1: Full-model table with stars appended to estimates
  full_tbl <- td_full %>%
    filter(term != "(Intercept)") %>%
    transmute(
      Term     = term,
      Estimate = paste0(round(estimate, 2), signif_stars(p.value))
    )
  # Append R² rows
  full_tbl <- rbind(
    full_tbl,
    data.frame(Term     = c("R-squared","Adj R-squared"),
               Estimate = c(as.character(r2), as.character(adjr2)),
               stringsAsFactors = FALSE)
  )
  full_grob <- tableGrob(full_tbl, rows = NULL, theme = tt_small)
  
  # Page 2: Simple-model table with stars on estimates
  simple_tbl <- lapply(predictor_names, function(pred) {
    m <- lm(as.formula(paste("y ~", pred)), data = df_tr)
    t <- tidy(m) %>% filter(term == pred)
    data.frame(Predictor = pred,
               Estimate   = paste0(round(t$estimate, 2), signif_stars(t$p.value)),
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
  simple_grob <- tableGrob(simple_tbl, rows = NULL, theme = tt_small)
  
  # Compute adjusted R² for each simple model
  simple_adj <- sapply(predictor_names, function(pred) {
    m <- lm(as.formula(paste("y ~", pred)), data = df_tr)
    round(summary(m)$adj.r.squared, 2)
  })
  simple_adj_tbl <- data.frame(
    Predictor = predictor_names,
    `Adj R-squared` = simple_adj,
    stringsAsFactors = FALSE
  )
  simple_adj_grob <- tableGrob(simple_adj_tbl, rows = NULL, theme = tt_small)
  
 
  preds <- list(
    Full       = predict(full_mod, df_te),
    Range      = predict(mod_range, df_te),
    Mean_Upper = predict(mod_mu, df_te),
    Ridge      = as.vector(model.matrix(y ~ ., df_te) %*% coef(mod_ridge))
  )
  p_pred <- imap(preds, ~ plot_pred_obs(df_te$y, .x, paste(.y, "–", score)))
  pred_grid <- do.call(grid.arrange, c(p_pred, ncol = 2))
  
  coef_plots <- list(Full = td_full, Range = td_range,
                     Mean_Upper = td_mu, Ridge = td_ridge)
  p_coef <- imap(coef_plots,
                 ~ plot_coefs(.x, paste(.y, "Coefficients –", score)))
  coef_grid <- do.call(grid.arrange, c(p_coef, ncol = 2))
  
  eff_plots <- lapply(keep, function(pred) {
    visreg(full_mod, xvar = pred, gg = TRUE,
           main = paste("Effect of", pred, "on", score))
  })
  effect_grid <- do.call(grid.arrange, c(eff_plots, ncol = 3))
  
  best2 <- td_full %>% filter(term != "(Intercept)") %>% arrange(p.value) %>% slice(1:2)
  best_plots <- lapply(best2$term, function(pred) {
    m   <- lm(as.formula(paste("y ~", pred)), data = df_tr)
    dfp <- df_te[, c("y", pred)]
    ggplot(dfp, aes_string(x = pred, y = "y")) +
      geom_point(size = 2) +
      geom_abline(slope = coef(m)[2], intercept = coef(m)[1], linetype = "dashed") +
      labs(x = pred, y = "Observed", title = paste("Simple Model: y ~", pred)) +
      theme_minimal()
  })
  best_grid <- do.call(grid.arrange, c(best_plots, ncol = 2))
  
  # Save per-score PDF
  out_pdf <- file.path(output_folder, paste0("full_summary_", score, ".pdf"))
  pdf(out_pdf, width = 14, height = 8.5)
  grid.arrange(full_grob,
               top = textGrob("Full Model: Estimates* (with significance stars)",
                              gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage()
  grid.arrange(simple_grob,
               top = textGrob("Simple Models: Estimates*", gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage()
  grid.arrange(simple_adj_grob,
               top = textGrob("Simple Models: Adjusted R-squared", gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage(); grid.arrange(pred_grid,
                               top = textGrob("Predicted vs Observed", gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage(); grid.arrange(coef_grid,
                               top = textGrob("Coefficient Plots", gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage(); grid.arrange(effect_grid,
                               top = textGrob("Effect Plots (Full Model)", gp = gpar(fontsize = 14, fontface = "bold")))
  grid.newpage(); grid.arrange(best_grid,
                               top = textGrob("Simple Models: Top 2 Predictors", gp = gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  message("Saved: ", out_pdf)
}

# Combined tables wide format
# Full-model estimates wide
full_wide <- map_df(full_models, function(mod) {
  td <- tidy(mod) %>% filter(term != "(Intercept)")
  data.frame(Term = td$term,
             Estimate = paste0(round(td$estimate, 2), signif_stars(td$p.value)),
             stringsAsFactors = FALSE)
}, .id = "Score") %>%
  pivot_wider(names_from = Score, values_from = Estimate)

# Simple-model estimates wide
simple_wide <- map_df(score_names, function(score) {
  df <- data.frame(y = combined_df[[score]], combined_df[predictor_names])
  map_df(predictor_names, function(pred) {
    td <- tidy(lm(as.formula(paste("y ~", pred)), data = df)) %>% filter(term == pred)
    data.frame(Predictor = pred,
               Estimate  = paste0(round(td$estimate, 2), signif_stars(td$p.value)),
               stringsAsFactors = FALSE)
  }) %>% mutate(Score = score)
}) %>% pivot_wider(names_from = Score, values_from = Estimate)

# R² wide: build matrix manually
r2_vals <- map_dbl(full_models, ~ round(summary(.x)$r.squared, 2))
adjr2_vals <- map_dbl(full_models, ~ round(summary(.x)$adj.r.squared, 2))
r2_mat <- rbind(r2 = r2_vals, adjr2 = adjr2_vals)
colnames(r2_mat) <- names(full_models)

# Simple-model adjusted R² wide
simple_adj_mat <- sapply(score_names, function(score) {
  df <- data.frame(y = combined_df[[score]], combined_df[predictor_names])
  sapply(predictor_names, function(pred) {
    round(summary(lm(as.formula(paste("y ~", pred)), data = df))$adj.r.squared, 2)
  })
})
rownames(simple_adj_mat) <- predictor_names

# Write combined wide tables to PDF
pdf(file.path(output_folder, "all_scores_wide.pdf"), width = 18, height = 8)
# Page 1: Full-model estimates
grid.arrange(tableGrob(full_wide, rows = NULL, theme = tt_small),
             top = textGrob("Full Estimates* with Stars", gp = gpar(fontsize=12, fontface="bold")))
# Page 2: Simple-model estimates
grid.newpage()
grid.arrange(tableGrob(simple_wide, rows = NULL, theme = tt_small),
             top = textGrob("Simple Estimates* with Stars", gp = gpar(fontsize=12, fontface="bold")))
# Page 3: R² & Adj R²
grid.newpage()
grid.arrange(tableGrob(r2_mat, rows = c("R-squared","Adj R-squared"), theme = tt_small),
             top = textGrob("R² & Adj R² (Full Models)", gp = gpar(fontsize=12, fontface="bold")))
# Page 4: Simple-model Adjusted R²
grid.newpage()
grid.arrange(tableGrob(simple_adj_mat, rows = predictor_names, theme = tt_small),
             top = textGrob("Adjusted R² for Simple Models", gp = gpar(fontsize=12, fontface="bold")))

simple_r2_mat <- sapply(score_names, function(score) {
  df <- data.frame(y = combined_df[[score]], combined_df[predictor_names])
  sapply(predictor_names, function(pred) {
    summary(lm(as.formula(paste("y ~", pred)), data = df))$r.squared
  })
})
rownames(simple_r2_mat) <- predictor_names

# compute average R² per predictor
avg_r2 <- rowMeans(simple_r2_mat)
avg_r2_df <- data.frame(
  Predictor = predictor_names,
  Avg_R2    = round(avg_r2, 2),
  stringsAsFactors = FALSE
)
# order factor levels so highest appears rightmost
avg_r2_df$Predictor <- factor(
  avg_r2_df$Predictor,
  levels = avg_r2_df$Predictor[order(avg_r2_df$Avg_R2)]
)

# build the bar‐plot with vertical bars and y‐limits 0 to 1
p_avg_r2 <- ggplot(avg_r2_df, aes(x = Predictor, y = Avg_R2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x     = NULL,
    y     = expression("Average " ~ R^2),
    title = "Average R² Across Scores for Simple Models"
  ) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# add new page and draw
grid.newpage()
grid.arrange(
  p_avg_r2,
  top = textGrob(
    "Average R² for Simple (Single-Predictor) Models",
    gp  = gpar(fontsize = 12, fontface = "bold")
  )
)
# --- NEW Page 6: Hierarchical Clustering of Scores ---
grid.newpage()

dist_scores <- dist(t(simple_r2_mat), method = "euclidean")

# perform clustering (you can swap in "average", "complete", or "ward.D2")
hc_scores <- hclust(dist_scores, method = "ward.D2")

# number of clusters
k <- 4

# get a distinct palette
cluster_cols <- RColorBrewer::brewer.pal(k, "Set1")

# plot the dendrogram
plot(
  hc_scores,
  main = "Score Clustering by Simple-Model R² Profiles",
  xlab = "",
  sub  = ""
)

# draw one box per cluster, in its own color
rect.hclust(
  hc_scores,
  k      = k,
  border = cluster_cols
)

# legend showing which color = which cluster
legend(
  "topright",
  legend = paste("Cluster", 1:k),
  fill   = cluster_cols,
  bty    = "n",
  cex    = 0.8
)

clusters   <- cutree(hc_scores, k = k)
sil        <- silhouette(clusters, dist_scores)

# Build a data.frame of per‐score silhouette widths
sil_df <- data.frame(
  Score        = names(clusters),
  Cluster      = clusters,
  Silhouette   = round(sil[, "sil_width"], 2),
  stringsAsFactors = FALSE
)
# Compute average silhouette by cluster
avg_sil_by_clust <- sil_df %>%
  group_by(Cluster) %>%
  summarise(Avg_Sil = round(mean(Silhouette), 2)) %>%
  arrange(Cluster)

# And overall average
overall_avg <- round(mean(sil_df$Silhouette), 2)

# Now add to your PDF device (after the dendrogram page):
grid.newpage()
grid.arrange(
  tableGrob(sil_df, rows = NULL, theme = tt_small),
  top = textGrob(
    paste0("Silhouette Widths (k = ", k, ")"),
    gp  = gpar(fontsize = 12, fontface = "bold")
  )
)
# And a small summary table beneath it
grid.newpage()
grid.arrange(
  tableGrob(
    rbind(
      as.data.frame(avg_sil_by_clust),
      data.frame(Cluster = "Overall", Avg_Sil = overall_avg)
    ),
    rows = NULL,
    theme = tt_small
  ),
  top = textGrob(
    "Average Silhouette Width by Cluster",
    gp  = gpar(fontsize = 12, fontface = "bold")
  )
)
# --- NEW: permutation tests (overall & per-cluster) for silhouette means ---
set.seed(123)
n_perm <- 1000

# observed per-cluster and overall means (you already have these)
obs_cluster_means <- avg_sil_by_clust$Avg_Sil
obs_overall_mean  <- overall_avg  # from earlier

# run permutations: shuffle labels, compute cluster means
perm_cluster_mat <- replicate(n_perm, {
  labs_perm <- sample(clusters)
  sil_p     <- silhouette(labs_perm, dist_scores)[, "sil_width"]
  tapply(sil_p, labs_perm, mean)
})
# perm_cluster_mat is a k × n_perm matrix

# per-cluster p-values
p_cluster <- sapply(seq_len(nrow(perm_cluster_mat)), function(i) {
  mean(perm_cluster_mat[i, ] >= obs_cluster_means[i])
})

# overall silhouette mean under each perm: rowMeans of perm_cluster_mat
perm_overall_means <- colMeans(perm_cluster_mat)
p_overall <- mean(perm_overall_means >= obs_overall_mean)

# assemble into one table
perm_summary <- data.frame(
  Cluster        = c("Overall", paste0("Cluster ", avg_sil_by_clust$Cluster)),
  Observed_Mean  = round(c(obs_overall_mean, obs_cluster_means), 2),
  Permutation_p  = round(c(p_overall, p_cluster),      3),
  stringsAsFactors = FALSE
)

# add a new page with this summary
grid.newpage()
grid.arrange(
  tableGrob(perm_summary, rows = NULL, theme = tt_small),
  top = textGrob(
    paste0("Permutation Test for Mean Silhouette (n = ", n_perm, ")"),
    gp  = gpar(fontsize = 12, fontface = "bold")
  )
)

# now close the PDF
dev.off()
message("Saved all_scores_wide.pdf with permutation p-value for clustering") 