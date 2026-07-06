# =============================================================================
# correlation_resolution.R — correlation vs sampling resolution
# =============================================================================
# WHAT IT DOES: Plots how index correlations change with the sampling resolution of the
#   reaction norm.
# REQUIRES:     Regression-summary CSVs under output/regression_summary_stats/.
# PRODUCES:     Correlation-vs-resolution figure(s).
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","correlation_resolution.R"))
# -----------------------------------------------------------------------------
# PARAMETERS (edit the named assignment near the top of the script)
#   directory    path, default output/regression_summary_stats  input CSV folder          [COMMON]
#   resolutions  integer vector, default c(50,25,20,15,10,5,2,1) sampling resolutions
#   n_samples    integer vector, default c(2,3,4,5,6,11,26,50)   sample counts
#   sel_idx      integer vector, default c(1,3,7,8)   which resolutions to plot            [COMMON]
#   files        character vector  the input CSVs
#   stat_cols    character vector  summary-stat columns to use
# =============================================================================

library(tidyverse)
library(patchwork)
library(viridis)
library(lme4)
library(lmerTest)
library(ggplot2)
library(forcats)
library(dplyr)   # must come last — re-asserts dplyr::select/count over MASS

directory   <- here::here("output", "regression_summary_stats")
resolutions <- c(50, 25, 20, 15, 10, 5, 2, 1)
n_samples   <- c( 2,  3, 4,  5,   6, 11,26,50)
files       <- c(
  "regression_data_full_interval_49_indices_1_50.csv",
  "regression_data_full_interval_25_indices_1_26_50.csv",
  "regression_data_full_interval_20_indices_1_21_41_50.csv",
  "regression_data_full_interval_15_indices_1_16_31_46_50.csv",
  "regression_data_full_interval_10_indices_1_11_21_31_41_50.csv",
  "regression_data_full_interval_5_indices_1_6_11_16_21_26_31_36_41_46_50.csv",
  "regression_data_full_interval_2_indices_1_3_5_7_9_11_13_15_17_19_21_23_25_27_29_31_33_35_37_39_41_43_45_47_49_50.csv",
  "regression_data_full_interval_1_indices_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34_35_36_37_38_39_40_41_42_43_44_45_46_47_48_49_50.csv"
)


sel_idx   <- c(1,3,7,8)
sel_files <- files[sel_idx]
sel_n     <- n_samples[sel_idx]


stat_cols <- c("min","max","mean","median","slope","range","variance","mean_lower","mean_upper")
df_all <- map2_df(sel_files, sel_n, ~{
  dat <- read_csv(file.path(directory, .x), show_col_types = FALSE)
  scores <- dat %>%
    dplyr::select(all_of(setdiff(colnames(dat), c(colnames(dat)[1], stat_cols)))) %>%
    pivot_longer(everything(), names_to="metric", values_to="value") %>%
    mutate(type="score", sample_size=.y)
  sums <- dat %>%
    dplyr::select(all_of(stat_cols)) %>%
    pivot_longer(everything(), names_to="metric", values_to="value") %>%
    mutate(type="summary", sample_size=.y)
  bind_rows(scores, sums)
})


score_metrics   <- df_all %>% 
  filter(type == "score") %>% 
  pull(metric) %>% 
  unique()

summary_metrics <- df_all %>% 
  filter(type == "summary") %>% 
  pull(metric) %>% 
  unique()

# drop any that also appear in score_metrics
summary_metrics <- setdiff(summary_metrics, score_metrics)

# now re-level so all scores come first, then the non-overlapping summaries
df_all <- df_all %>%
  mutate(
    metric = factor(metric,
                    levels = c(score_metrics, summary_metrics))
  )
# global y-limits
y_limits <- range(df_all$value, na.rm = TRUE)


combined=ggplot(df_all, aes(x = metric, y = value, fill = type)) +
  geom_boxplot(
    outlier.size = 0.7,
    color        = "black",
    position     = position_dodge(width = 0.8),
    width        = 0.6
  ) +
  scale_fill_manual(
    values = c(score = "steelblue", summary = "tomato"),
    name   = NULL
  ) +
  scale_y_continuous(limits = y_limits) +
  facet_wrap(
    ~ sample_size,
    ncol   = 2,
    scales = "free_x",
    labeller = labeller(sample_size = function(x) paste0("n = ", x))
  ) +
  labs(
    x     = NULL,
    y     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background   = element_rect(fill = "grey90"),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

path <- here::here("output", "plots")
if (!dir.exists(dirname(path))) {
  dir.create(dirname(path), recursive = TRUE)
}
ggsave(create.dir = TRUE, 
  combined,
  filename = file.path(here::here("output", "plots", "boxplot_resolution.pdf")),
  width    = 6.3,
  height   = 8,
  dpi      = 900
)

print(combined)


###############################################################################################################################################################################################

wide_list <- vector("list", length(files))
for(i in seq_along(files)){
  df <- read.csv(file.path(directory, files[i]), stringsAsFactors = FALSE)
  PS      <- df[, !colnames(df) %in% c(colnames(df)[1], stat_cols), drop = FALSE]
  summary <- df[, stat_cols, drop = FALSE]
  
  cors <- matrix(NA, nrow=ncol(PS)*ncol(summary), ncol=3,
                 dimnames=list(NULL, c("PS","Stat","r")))
  idx <- 1
  for(p in colnames(PS)){
    for(s in colnames(summary)){
      cors[idx, "PS"]   <- p
      cors[idx, "Stat"] <- s
      cors[idx, "r"]    <- cor(PS[[p]], summary[[s]], use="pairwise.complete.obs",method = "pearson")
      idx <- idx + 1
    }
  }
  wide_list[[i]] <- data.frame(cors, stringsAsFactors=FALSE)
  wide_list[[i]]$r <- as.numeric(wide_list[[i]]$r)
  names(wide_list)[i] <- paste0("n", n_samples[i])
}

wide_df <- wide_list[[1]]
wide_df <- setNames(wide_df, c("PS","Stat","n2"))

for(i in 2:length(wide_list)){
  df_i <- wide_list[[i]]
  wide_df <- merge(wide_df, df_i,
                   by=c("PS","Stat"), all=TRUE, suffixes=c("",paste0(".n",i)))
  
  colnames(wide_df)[ncol(wide_df)] <- paste0("n", n_samples[i])
}

nr        <- nrow(wide_df)
model_df  <- data.frame(
  PS        = wide_df$PS,
  Stat      = wide_df$Stat,
  intercept = NA_real_,
  slope     = NA_real_,
  se_slope  = NA_real_,
  nobs      = NA_integer_,
  df        = NA_integer_,
  t_stat    = NA_real_,
  p_value   = NA_real_,
  stringsAsFactors = FALSE
)


for(j in seq_len(nr)) {
  y   <- as.numeric(wide_df[j, paste0("n", n_samples)])
  x   <- n_samples
  ok  <- !is.na(y)
  n_j <- sum(ok)
  if(n_j >= 3) {              # need at least 3 points to get df≥1
    fit            <- lm(y[ok] ~ x[ok])
    cf             <- summary(fit)$coefficients
    model_df$intercept[j] <- cf["(Intercept)","Estimate"]
    model_df$slope[j]     <- cf["x[ok]","Estimate"]
    model_df$se_slope[j]  <- cf["x[ok]","Std. Error"]
    model_df$nobs[j]      <- n_j
    model_df$df[j]        <- n_j - 2
    # t‐stat and p‐value
    tval                  <- model_df$slope[j] / model_df$se_slope[j]
    model_df$t_stat[j]    <- tval
    model_df$p_value[j]   <- 2 * pt(-abs(tval), df = model_df$df[j])
  }
}


#################################### visualisation of correlation resolution ##########################################

delta <- 0.001  # Define the equivalence margin for slope values.
# This sets the threshold below which slopes are considered 
# equal to zero 

model_df$peq <- NA_real_

for (j in seq_len(nrow(model_df))) {
  b   <- model_df$slope[j]
  se  <- model_df$se_slope[j]
  df  <- model_df$df[j]
  if (!is.na(b) && !is.na(se) && df > 0) {
    # test 1: H0: beta <= -delta  vs  HA: beta > -delta
    t1 <- (b - (-delta)) / se
    
    p1 <- pt(t1, df, lower.tail = FALSE)
    
    # test 2: H0: beta >= +delta  vs  HA: beta < +delta
    t2 <- (delta - b) / se
    
    
    p2 <- pt(t2, df, lower.tail = FALSE)
  
    # equivalence p‐value is the maximum of the two one‐sided p’s
    model_df$peq[j] <- max(p1, p2)
  }
}

# now “equivalent to zero” at alpha<0.05 --> significance value of the equivalence test
model_df$equiv_robust <- model_df$peq < 0.05


robust_eq <- subset(model_df, equiv_robust)
print(robust_eq[, c("PS","Stat","slope","se_slope","df","peq")])

library(tidyverse)

###### classification dataframe and plot 
class_df <- model_df %>%
  transmute(
    PS, 
    Stat, 
    flag = if_else(equiv_robust, "robust", "sensitive")
  )


heat_df <- class_df %>%
  pivot_wider(
    names_from  = Stat,
    values_from = flag,
    values_fill = "sensitive"
  )



hm_long <- heat_df %>%
  pivot_longer(
    cols      = -PS,        # everything except PS
    names_to  = "Stat",
    values_to = "flag"
  ) %>%
  mutate(
    PS   = factor(PS,   levels = unique(PS)),
    Stat = factor(Stat, levels = unique(Stat))
  )


p_heatmap <- ggplot(hm_long, aes(x = Stat, y = PS, fill = flag)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(robust = "forestgreen", sensitive = "firebrick"),
    name   = "Classification"
  ) +
  labs(
    x     = "Summary Statistic",
    y     = "Plasticity Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 8),
    legend.position = "bottom"
  )
print(p_heatmap)


summary_counts <- hm_long %>%
  dplyr::count(PS, flag) %>%
  pivot_wider(names_from = flag, values_from = n, values_fill = 0) %>%
  mutate(
    total   = robust + sensitive,
    prop    = robust / total
  )

print(summary_counts)

p_summary <- ggplot(summary_counts, aes(x = reorder(PS, prop), y = prop)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Plasticity Score",
    y = "% robust"
  ) +
  theme_minimal(base_size = 10)

print(p_summary)

###############################################################
# Heatmap: slope of correlation ~ sample_size, labelled with p-values
###############################################################

slope_heat <- model_df %>%
  filter(!is.na(slope)) %>%
  mutate(
    p_label = ifelse(is.na(p_value), "",
                ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*", "")))),
    PS   = factor(PS,   levels = rev(sort(unique(PS)))),
    Stat = factor(Stat, levels = sort(unique(Stat)))
  )

p_slope_heat <- ggplot(slope_heat, aes(x = Stat, y = PS, fill = slope)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = p_label), size = 2.8, colour = "black") +
  scale_fill_gradient2(
    low      = "#4575B4",
    mid      = "white",
    high     = "#D73027",
    midpoint = 0,
    name     = "Slope"
  ) +
  labs(
    x = "Summary Statistic",
    y = "Plasticity Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 9),
    legend.position = "right",
    panel.grid      = element_blank()
  )

print(p_slope_heat)

ggsave(create.dir = TRUE, 
  filename = here::here("output", "plots", "slope_heatmap_resolution.pdf"),
  plot     = p_slope_heat,
  width    = 6.3,
  height   = 7,
  dpi      = 900,
  device   = "pdf"
)
