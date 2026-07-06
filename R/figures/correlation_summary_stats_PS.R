# =============================================================================
# correlation_summary_stats_PS.R — summary-stat vs plasticity-score correlations
# =============================================================================
# WHAT IT DOES: Renders correlation panels between the reaction-norm summary statistics
#   and the plasticity scores.
# REQUIRES:     A regression-data CSV under output/regression_summary_stats/.
# PRODUCES:     Correlation-panel figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","correlation_summary_stats_PS.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   the input CSV path (read via read_csv at the top)                                      [COMMON]
#   stat_cols   character vector  summary-stat columns to use
# =============================================================================

# install.packages(c(
#   "readr","dplyr","reshape2","ggplot2","scales","Hmisc","broom","RColorBrewer"
# ))
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(Hmisc)      # rcorr()
library(broom)      # tidy()
library(RColorBrewer)
library(psych)
library(dplyr)
library(igraph)
library(ggraph)
library(tidyr)
library(ggforce)
library(ggalluvial)
library(tidyverse)
library(patchwork)
library(viridis)
library(factoextra)


############################# data reading and computation of the linear regression for the slope heatmaps 
df <- read_csv(here::here("data", "regression_summary_stats", "regression_data_full_interval_10_indices_1_11_21_31_41_50.csv"))


stat_cols <- c("min","max","mean","median","slope","range","variance","mean_lower","mean_upper")
PS      <- df[, !colnames(df) %in% c(colnames(df)[1], stat_cols), drop = FALSE]
sumstat <- df[, stat_cols, drop = FALSE]

ct <- corr.test(
  as.matrix(cbind(PS, sumstat)),
  method = "pearson",
  adjust = "none"         # turn off p-value adjustment
)

r_mat <- ct$r[colnames(PS), colnames(sumstat)]
p_mat <- ct$p[colnames(PS), colnames(sumstat)]
PS_names   <- colnames(PS)
Stat_names <- colnames(sumstat)


slope_mat     <- matrix(NA,
                        nrow = length(PS_names),
                        ncol = length(Stat_names),
                        dimnames = list(PS_names, Stat_names))

cor_df <- melt(r_mat, varnames = c("PS","Stat"), value.name = "r") %>%
  left_join(
    melt(p_mat, varnames = c("PS","Stat"), value.name = "p"),
    by = c("PS","Stat")
  ) %>%
  mutate(star = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  ))

pval_mat      <- slope_mat   # same dims for the p‐values
slope_ci_lo   <- slope_mat   # lower CIs
slope_ci_hi   <- slope_mat   # upper CIs

for(ps in PS_names) {
  for(st in Stat_names) {
    sub <- na.omit(df[, c(ps, st)])
    fm  <- lm(sub[[ps]] ~ sub[[st]])
    td  <- broom::tidy(fm, conf.int = TRUE)
    
    slope_mat[ps, st]     <- td$estimate[2]
    pval_mat[ps, st]      <- td$p.value[2]
    slope_ci_lo[ps, st]   <- td$conf.low[2]
    slope_ci_hi[ps, st]   <- td$conf.high[2]
  }
}


slope_ci_w  <- slope_ci_hi - slope_ci_lo
max_w       <- max(slope_ci_w, na.rm = TRUE)
precision_mat <- 1 - (slope_ci_w / max_w)


slope_df <- reshape2::melt(slope_mat,      varnames = c("PS","Stat"), value.name = "slope") %>%
  left_join(
    reshape2::melt(pval_mat,      varnames = c("PS","Stat"), value.name = "pval"),
    by = c("PS","Stat")
  ) %>%
  left_join(
    reshape2::melt(precision_mat, varnames = c("PS","Stat"), value.name = "precision"),
    by = c("PS","Stat")
  ) %>%
  mutate(
    star = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )


####################################################### just the theme for some of the figures 
theme_prof <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title       = element_text(face = "bold"),
    legend.title     = element_text(face = "bold"),
    legend.key.size  = unit(0.8, "lines")
  )

########################################################## heatmap of the correlation between the summary stats and the scores --> i really liked this, Zoran doesnt 
# 7) Heatmap A: Correlation
p_corr <- ggplot(cor_df, aes(x = Stat, y = PS, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = brewer.pal(9, "Blues")[5],
    mid      = "white",
    high     = brewer.pal(9, "Reds")[5],
    midpoint = 0,
    limits   = c(-1,1),
    name     = expression(r)
  ) +
  geom_text(aes(label = sprintf("%.2f%s", r, star)), size = 3) +
  labs(
    x = "Summary Statistic",
    y = "Plasticity Score"
  ) +
  theme_prof


print(p_corr)
#ggsave(create.dir = TRUE, "heatmap_correlation.png", p_corr, width = 12, height = 9, dpi = 300)


######################################################## Heatmap containing the slopes of the single predictor regression model (one model per score ~ summary stat combination) 
#--> later there is also a multi predictor regression model. The opacity of the tiles indicates the confidence interval width, the stars indicate significance
# Zoran ruled this out bc it looks too crowded and additinally linear reghression with one preductor is simply correlation but with offset. This can therefor be left out 

p_slope <- ggplot(slope_df, aes(x = Stat, y = PS, fill = precision)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low  = "white",
    high = brewer.pal(9, "YlGnBu")[5],
    name = "Precision"
  ) +
  geom_text(aes(label = sprintf("%.2f%s", slope, star)),
            size = 3, color = "black") +
  labs(
    x = "Summary Statistic",
    y = "Plasticity Score"
  ) +
  theme_prof +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(p_slope)
#ggsave(create.dir = TRUE, here::here("output", "plots", "scatter_plots", "heatmap_slope_with_stars.png"), p_slope, width = 12, height = 9, dpi = 300)

################################################################ this creates scatter plots for each and every combination of summary stat and plasticity score --> linear regression plot 
# I didnt use this anywhere but it is useful for looking the regressions up 
# another interesting thing might be that there are rather interesting patterns in the data sometimes with which one can distinguish between different genotype forms (linear genotype form correlates with the variance differently than gaussian for example and this is really nicely visible in those plots)

# 9) Scatter + LM + annotations: helper function
plot_scatter <- function(data, xvar, yvar) {
  sub <- na.omit(data[, c(xvar, yvar)])
  fm  <- lm(sub[[yvar]] ~ sub[[xvar]])
  sm  <- summary(fm)
  r   <- cor(sub[[xvar]], sub[[yvar]])
  r2  <- sm$r.squared
  b1  <- coef(fm)[2]
  ci  <- confint(fm)[2, ]
  p   <- coef(sm)[2, "Pr(>|t|)"]
  
  ggplot(sub, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.8) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = sprintf(
        "r = %.2f\nR² = %.2f\nβ₁ = %.2f (%.2f, %.2f)\np = %.2g",
        r, r2, b1, ci[1], ci[2], p
      ),
      hjust = 1.05, vjust = -0.3,
      size = 3.2
    ) +
    labs(
      x = xvar,
      y = yvar
    ) +
    theme_prof
}

# 10) Generate & save scatter plots for every PS × summary-stat pair
out_base <- here::here("output", "plots", "scatter_plots")
if (!dir.exists(out_base)) {
  dir.create(out_base, recursive = TRUE)
}
for(ps in colnames(PS)) {
  dir.create(file.path(out_base, ps), recursive = TRUE, showWarnings = FALSE)
  for(st in colnames(sumstat)) {
    p_sc <- plot_scatter(df, st, ps)
    ggsave(create.dir = TRUE, 
      filename = file.path(out_base, ps, paste0(ps, "_vs_", st, ".png")),
      plot     = p_sc,
      width    = 6, height = 5, dpi = 300
    )
  }
}

##########################################################

make_dend_plot <- function(hc, title, k = 3, shrink = 0.5) {
  # rebuild the shrink transform
  shrink_trans <- scales::trans_new(
    name      = paste0("shrink_", shrink),
    transform = function(x) x * shrink,
    inverse   = function(x) x / shrink
  )
  
  # fviz_dend with coloured branches
  p <- fviz_dend(hc,
                 k           = k,
                 rect        = FALSE,
                 show_labels = TRUE,
                 repel       = FALSE,
                 k_colors    = "jco",
                 cex         = 0.8
  ) +
    scale_y_continuous(trans = shrink_trans) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      plot.margin  = margin(t = 5, r = 5, b = 15, l = 5)
    )
  return(p)
}

#  HIERARCHICAL CLUSTERING OF CORRELATIONS 

d_ps <- dist(r_mat)
hc_ps <- hclust(d_ps, method = "complete")
# Plot with coloured branches, k = 3 clusters
p_dend_ps <- make_dend_plot(hc_ps, "Clustering of plasticity scores with regards to correlation withsummary stat ", k = 3)
print(p_dend_ps)
# one can yet again confirm that the clustering is in agreement with the other previous clusterings but there might also be no point in that 


################################################## multiple predictor regression heatmap --> range and mean_upper yield NA since there are colinearities for those

n_ps <- ncol(PS)
n_st <- ncol(sumstat)
multi_coef_full  <- matrix(NA, n_ps, n_st, dimnames = list(colnames(PS), colnames(sumstat)))
multi_ci_lo_full <- multi_coef_full
multi_ci_hi_full <- multi_coef_full
multi_pval_full  <- multi_coef_full

# Loop over each Plasticity Score
for(ps in colnames(PS)) {
  dat  <- na.omit(cbind(PS = PS[[ps]], sumstat))
  fit  <- lm(PS ~ ., data = dat)
  td   <- broom::tidy(fit, conf.int = TRUE) %>% filter(term != "(Intercept)")
  for(i in seq_len(nrow(td))) {
    st <- td$term[i]
    multi_coef_full[ps, st]  <- td$estimate[i]
    multi_ci_lo_full[ps, st] <- td$conf.low[i]
    multi_ci_hi_full[ps, st] <- td$conf.high[i]
    multi_pval_full[ps, st]  <- td$p.value[i]
  }
}

# Compute precision for full model
multi_ci_w_full <- multi_ci_hi_full - multi_ci_lo_full
max_w_full     <- max(multi_ci_w_full, na.rm = TRUE)
multi_prec_full <- 1 - (multi_ci_w_full / max_w_full)

# Melt into data frame
multi_df_full <- reshape2::melt(multi_coef_full,   varnames = c("PS","Stat"), value.name = "beta") %>%
  left_join(reshape2::melt(multi_prec_full, varnames = c("PS","Stat"), value.name = "precision"), by = c("PS","Stat")) %>%
  left_join(reshape2::melt(multi_pval_full,  varnames = c("PS","Stat"), value.name = "pval"), by = c("PS","Stat")) %>%
  mutate(star = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ ""
  ))

# Plot full heatmap: β fill, alpha ∝ precision
p_multi_full <- ggplot(multi_df_full, aes(x = Stat, y = PS)) +
  geom_tile(aes(fill = beta, alpha = precision), color = "white") +
  scale_fill_gradient2(low = RColorBrewer::brewer.pal(9, "RdBu")[1],
                       mid = "white",
                       high = RColorBrewer::brewer.pal(9, "RdBu")[9],
                       midpoint = 0,
                       name = expression(beta)) +
  scale_alpha_continuous(range = c(0.3,1), name = "Precision") +
  geom_text(aes(label = sprintf("%.2f%s", beta, star)), size = 3) +
  labs(x = "Summary Statistic", y = "Plasticity Score") +
  theme_prof +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_multi_full)

########################################

# --- R² FOR FULL MULTI-PREDICTOR MODELS ----
r2_full <- sapply(colnames(PS), function(ps) {
  dat <- na.omit(cbind(PS = PS[[ps]], sumstat))
  summary(lm(PS ~ ., data = dat))$r.squared
})
r2_full_df <- data.frame(PS = names(r2_full), R2 = r2_full)
p_r2_full <- ggplot(r2_full_df, aes(x = PS, y = R2)) +
  geom_col(fill = "steelblue") +
  labs(x = "Plasticity Score", y = expression(R^2)) +
  theme_prof +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_r2_full)


############################################### multiple predictor regression heatmap with the collinearities removed such that we dont get NA for range and mean_upper

substats_sub <- sumstat[, 3:8] 
n_st_sub    <- ncol(substats_sub)
multi_coef_sub  <- matrix(NA, n_ps, n_st_sub, dimnames = list(colnames(PS), colnames(substats_sub)))
multi_ci_lo_sub <- multi_coef_sub
multi_ci_hi_sub <- multi_coef_sub
multi_pval_sub  <- multi_coef_sub

for(ps in colnames(PS)) {
  dat <- na.omit(cbind(PS = PS[[ps]], substats_sub))
  fit <- lm(PS ~ ., data = dat)
  td  <- broom::tidy(fit, conf.int = TRUE) %>% filter(term != "(Intercept)")
  for(i in seq_len(nrow(td))) {
    st <- td$term[i]
    multi_coef_sub[ps, st]  <- td$estimate[i]
    multi_ci_lo_sub[ps, st] <- td$conf.low[i]
    multi_ci_hi_sub[ps, st] <- td$conf.high[i]
    multi_pval_sub[ps, st]  <- td$p.value[i]
  }
}

multi_ci_w_sub <- multi_ci_hi_sub - multi_ci_lo_sub
max_w_sub     <- max(multi_ci_w_sub, na.rm = TRUE)
multi_prec_sub <- 1 - (multi_ci_w_sub / max_w_sub)

multi_df_sub <- reshape2::melt(multi_coef_sub, varnames = c("PS","Stat"), value.name = "beta") %>%
  left_join(reshape2::melt(multi_prec_sub, varnames = c("PS","Stat"), value.name = "precision"), by = c("PS","Stat")) %>%
  left_join(reshape2::melt(multi_pval_sub, varnames = c("PS","Stat"), value.name = "pval"), by = c("PS","Stat")) %>%
  mutate(star = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ ""
  ))

# Plot subset heatmap: β fill, alpha ∝ precision
p_multi_sub <- ggplot(multi_df_sub, aes(x = Stat, y = PS)) +
  geom_tile(aes(fill = beta, alpha = precision), color = "white") +
  scale_fill_gradient2(low = RColorBrewer::brewer.pal(9, "YlGnBu")[1],
                       mid = "white",
                       high = RColorBrewer::brewer.pal(9, "YlGnBu")[5],
                       midpoint = 0,
                       name = expression(beta)) +
  scale_alpha_continuous(range = c(0.3,1), name = "Precision") +
  geom_text(aes(label = sprintf("%.2f%s", beta, star)), size = 3) +
  labs(x = "Summary Statistic", y = "Plasticity Score") +
  theme_prof +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_multi_sub)

################################# R² FOR SUBSET MULTI-PREDICTOR MODELS 

r2_sub <- sapply(colnames(PS), function(ps) {
  dat <- na.omit(cbind(PS = PS[[ps]], substats_sub))
  summary(lm(PS ~ ., data = dat))$r.squared
})
r2_sub_df <- data.frame(PS = names(r2_sub), R2 = r2_sub)
p_r2_sub <- ggplot(r2_sub_df, aes(x = PS, y = R2)) +
  geom_col(fill = "steelblue") +
  labs(x = "Plasticity Score", y = expression(R^2)) +
  theme_prof +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p_r2_sub)


###########################################
# --- PRECISION-ONLY HEATMAP FOR FULL MODEL ----
p_multi_precision <- ggplot(multi_df_full, aes(x = Stat, y = PS, fill = precision)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white",
                      high = brewer.pal(9, "YlGnBu")[5],
                      name = "Precision") +
  geom_text(aes(label = sprintf("%.2f%s", beta, star)), size = 3) +
  labs(x = "Summary Statistic", y = "Plasticity Score") +
  theme_prof
print(p_multi_precision)

# --- PRECISION-ONLY HEATMAP FOR SUBSET MODEL ----
p_multi_precision_sub <- ggplot(multi_df_sub, aes(x = Stat, y = PS, fill = precision)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white",
                      high = brewer.pal(9, "YlGnBu")[5],
                      name = "Precision") +
  geom_text(aes(label = sprintf("%.2f%s", beta, star)), size = 3) +
  labs(x = "Summary Statistic", y = "Plasticity Score") +
  theme_prof
print(p_multi_precision_sub)


######################################
final_fig3 <- wrap_plots(
  A = wrap_elements(full = p_corr),
  B = wrap_elements(full = p_slope),
  C= wrap_elements(full = p_dend_ps),
  design = "
AAABBB
AAABBB
AAABBB
CCCCCC

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


out= here::here("output", "plots", "figures")
if (!dir.exists(out)) {
  dir.create(out, recursive = TRUE)
}
ggsave(create.dir = TRUE, 
  filename = here::here("output", "plots", "figures", "sum_stats_cor_pearson.pdf"),
  plot     = final_fig3,
  device   = "pdf",
  width    = 6.3,
  height   = 8,
  units    = "in",
  dpi      = 900
)

ggsave(create.dir = TRUE, 
  filename = here::here("output", "plots", "figures", "hclust_dend_corrdist.pdf"),
  plot     = p_dend_ps,
  device   = "pdf",
  width    = 6.3,
  height   = 5,
  units    = "in",
  dpi      = 900
)


######################################

final_fig4 <- wrap_plots(
  A = wrap_elements(full = p_multi_precision),
  B = wrap_elements(full = p_multi_precision_sub),
  C= wrap_elements(full = p_r2_full),
  D = wrap_elements(full = p_r2_sub),
  design = "
AAABBB
AAABBB
AAABBB
CCCDDD
CCCDDD
"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

print(final_fig4)

ggsave(create.dir = TRUE, 
  filename = here::here("output", "plots", "figures", "ols.pdf"),
  plot     = final_fig4,
  device   = "pdf",
  width    = 6.3,
  height   = 8,
  units    = "in",
  dpi      = 900
)



########################################## these are the figures after the corrections that zoran suggested

########################################## alluv plot that is in the figure

thr_r <- 0.55
thr_p <- 0.05


df <- cor_df %>%
  filter(abs(r) >= thr_r, p < thr_p) %>%
  mutate(
    id      = row_number(),
    abs_rho = abs(r)
  )

pos1 <- df %>% 
  filter(r > 0) %>%
  transmute(id, abs_rho,
            x       = "Positive",
            stratum = Stat,
            flow_stat = Stat)

pos2 <- df %>% 
  filter(r > 0) %>%
  transmute(id, abs_rho,
            x       = "Score",
            stratum = PS,
            flow_stat = Stat)

neg1 <- df %>% 
  filter(r < 0) %>%
  transmute(id, abs_rho,
            x       = "Score",
            stratum = PS,
            flow_stat = Stat)

neg2 <- df %>% 
  filter(r < 0) %>%
  transmute(id, abs_rho,
            x       = "Negative",
            stratum = Stat,
            flow_stat = Stat)

df_lodes <- bind_rows(pos1, pos2, neg1, neg2) %>%
  mutate(
    x = factor(x, levels = c("Positive","Score","Negative"))
  )


alluv=ggplot(df_lodes,
       aes(x        = x,
           alluvium = id,
           y        = abs_rho,
           stratum  = stratum,
           fill     = flow_stat,
           alpha    = abs_rho)) +
  # ribbons
  geom_flow(stat   = "flow",
            width  = 0.4,
            colour = "grey20") +
  # left/right stat bars
  geom_stratum(data   = df_lodes %>% filter(x != "Score"),
               aes(fill = stratum),
               width  = 0.4,
               colour = "grey30") +
  # middle PS bars
  geom_stratum(data  = df_lodes %>% filter(x == "Score"),
               aes(stratum = stratum),
               fill   = "grey90",
               colour = "grey30",
               width  = 0.4) +
  # labels in every bar
  geom_text(stat   = "stratum",
            aes(label = after_stat(stratum)),
            size   = 4,
            colour = "black") +
  # remove the fill‐(Stat) legend
  scale_fill_viridis_d(option = "plasma", 
                       begin  = 0.4,    # skip the darkest 15%
                       end    = 0.9,    # skip the very brightest 15%
                       guide  = "none") +
  # show only an alpha‐legend for |r|
  scale_alpha_continuous(range = c(0.3, 1),
                         name  = "Correlation strength\n(|r|)",
                         guide = guide_legend()) +
  scale_x_discrete(
    limits = c("Positive","Score","Negative"),
    labels = c("r > 0.5", "Plasticity Scores", "r < -0.5"),
    expand = expansion(add = c(-1,-1))    # <-- adds one unit of padding on left and right
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "bottom"
  )


################################### "heatmap" with filtered correlations that are smaller or bigger than r = |0.55| --> not used right now 


ggplot(df, aes(x = Stat, y = PS)) +
  
  geom_point(aes(fill  = r,
                 alpha = abs_rho),
             shape  = 21,
             size   = 8,
             colour = "grey30") +
  # p-value stars on top
  geom_text(aes(label = star),
            colour = "black",
            size   = 3) +
  # diverging colour scale for sign & magnitude of r
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0,
    name     = "r"
  ) +
  # opacity encodes |r|
  scale_alpha_continuous(
    range = c(0.3, 1),
    name  = expression(r)
  ) +
  labs(
    title = "Filtered Correlations (|r| ≥ 0.55, p < 0.05)",
    x     = "Summary Statistic",
    y     = "Plasticity Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )





##################################################### this is the heatmap of the slopes of the single predictor regression model (one model per score ~ summary stat combination) (the reworked version that looks less crowded than the one on the top of this script )





# 1. Keep only significant rows (at least one star)
slope_df_sig <- slope_df %>%
  filter(star != "")

# 2. Plot with size = precision, fill = slope, and stars in the center
woho=ggplot(slope_df_sig, aes(x = Stat, y = PS)) +
  geom_point(aes(fill  = slope,
                 size  = precision),
             shape  = 21,        # circle with fill + border
             colour = "grey30",
             alpha  = 1) +     # fixed alpha for softness
  geom_text(aes(label = star),
            colour = "black",
            size   = 4) +
  scale_fill_gradient2(
    low      = "steelblue",   # negative slopes
    mid      = "white",       # around zero
    high     ="firebrick",  
    midpoint = 0,
    name     = expression(beta[1])
  ) +
  scale_size_continuous(
    range = c(5, 12),         # adjust min/max point size
    name  = "Precision"
  ) +

  
  labs(
    x = "Summary Statistic",
    y = "Plasticity Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.box.just = "center",
    plot.margin     = margin(t = 10, r = 5, b = 10, l = 5, unit = "pt")
  )

ggsave(create.dir = TRUE, here::here("output", "plots", "figures", "single_predictor_regression.pdf"), woho, width = 6.3, height = 8, units = "in", dpi = 900)



#
#final_fig5 <- wrap_plots(
#  A = wrap_elements(full = alluv),
#  B = wrap_elements(full = woho),
#  design = "
#AAABBBB
#AAABBBB
#AAABBBB
#AAABBBB
#"
#) +
#  plot_annotation(
#    tag_levels = "A",
#    theme      = theme(
#      plot.title        = element_text(size = 12, face = "bold"),
#      plot.tag.position = c(1, 1)
#    )
#  )
#
#print(final_fig5)
#
#ggsave(create.dir = TRUE, 
#  filename = here::here("output", "plots", "figures", "fig4.pdf"),
#  plot     = final_fig5,
#  device   = "pdf",
#  width    = 15,
#  height   = 15,
#  units    = "in",
#  dpi      = 900
#)

####################################### this is for the figure containing the hetmaps for the regression slopes for the full/subset model and not the single predictor model. Only the significant slopes are being reported. 
# withthis we can potentially make the claim that when using the summary statsitics in combination we can reproduce the scores by using the respecitve coefficients. In the end it is just proving that there certainly is a linear relationship between all the summary statisitcs and the plasticity scores (if one didnt believe it until now ;))
library(cowplot)   # for plot_grid
library(scales) 

df_full <- multi_df_full %>% mutate(Model = "Full")
df_sub  <- multi_df_sub  %>% mutate(Model = "Subset")

# filter to significant rows
sig_full <- multi_df_full %>% filter(star != "")
sig_sub  <- multi_df_sub  %>% filter(star != "")

make_dot <- function(df, title) {
  ggplot(df, aes(x = Stat, y = PS)) +
    geom_point(aes(fill = beta, size = precision),
               shape = 21, colour = "grey30") +
    geom_text(aes(label = star), size = 3) +
    scale_fill_gradient2(
      low      = brewer.pal(9, "RdBu")[1],
      mid      = "white",
      high     = brewer.pal(9, "RdBu")[9],
      midpoint = 0,
      name     = expression(beta)
    ) +
    scale_size_continuous(
      range   = c(2, 8),
      name    = "Precision",
      limits  = c(0, 1),
      breaks  = c(0, 0.5, 1),
      labels  = c("0.00", "0.50", "1.00")
    ) +
    guides(
      fill = guide_colourbar(order = 1, title.position = "top"),
      size = guide_legend(order = 2, title.position = "top")
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.box      = "vertical"
    )
}
sig_full <- multi_df_full %>% filter(star != "")
sig_sub  <- multi_df_sub  %>% filter(star != "")

p_full_sig <- make_dot(sig_full, "Full Model regression slopes")
p_sub_sig  <- make_dot(sig_sub,  "Subset Model regression slopes")

#––– 2) Mirrored R² plot with centered title + margins –––––––––––––––––––

df_r2 <- data.frame(
  PS        = names(r2_full),
  full_R2   = as.numeric(r2_full),
  subset_R2 = as.numeric(r2_sub),
  stringsAsFactors = FALSE
) %>%
  arrange(full_R2) %>%
  mutate(PS = factor(PS, levels = PS))

Rsqfull <- ggplot(df_r2) +
  geom_col(aes(x = -full_R2, y = PS),   fill = "steelblue") +
  geom_col(aes(x =  subset_R2, y = PS), fill = "firebrick") +
  geom_text(aes(x = 0, y = PS, label = PS),
            hjust = 0.5, size = 3, color = "black") +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", abs(x)),
    breaks = seq(-1, 1, by = 0.2),
    limits = c(-1, 1)
  ) +
  labs(
    x = expression(R^2),
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(5, 80, 5, 80)
  )

#––– 3) Combine into final figure –––––––––––––––––––––––––––––––––––––––––

final_fig6 <- wrap_plots(
  A = wrap_elements(full = p_full_sig),
  B = wrap_elements(full = p_sub_sig),
  C = wrap_elements(full = Rsqfull),
  design = "
AAABBB
AAABBB
AAABBB
CCCCCC
CCCCCC
"
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag          = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

print(final_fig6)
ggsave(create.dir = TRUE, 
  filename = here::here("output", "plots", "figures", "fig5.pdf"),
  plot     = final_fig6,
  device   = "pdf",
  width    = 6.3,
  height   = 8,
  units    = "in",
  dpi      = 900
)