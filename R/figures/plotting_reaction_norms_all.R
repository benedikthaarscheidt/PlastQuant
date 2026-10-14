# =============================================================================
# plotting_reaction_norms_all.R — reaction-norm overview figure
# =============================================================================
# WHAT IT DOES: Simulates a population with the model in R/sim (independent of the
#   heritability scenarios) and assembles the reaction-norm overview figure: per shape
#   the mother reaction norms of selected genotypes (smooth curve from the mother
#   parameters), their replicate measurements (dashed) and replicate range (ribbon),
#   plus histograms of the measured trait values.
# REQUIRES:     Nothing upstream (simulates its own data; no GWAS).
# PRODUCES:     output/plots/reaction_norms_final_figure.pdf,
#               reaction_norms_final_figure_panels.pdf, p_combined.pdf.
# HOW TO RUN:   setwd("~/PlastQuant"); source(here::here("R","figures","plotting_reaction_norms_all.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   H2              numeric, default 0.4    SNP share of the latent variance           [COMMON]
#   CAUSAL_SNP_NUM  integer, default 5      causal SNPs per parameter
#   NUM_GENOTYPES   integer, default 800    genotypes (divisible by 4, one block per shape)
#   NUM_SNPs        integer, default 2000   simulated SNPs (only the curves are used here)
#   SEED            integer, default 1
#   PLOT_GENOTYPES  list shape -> genotype IDs, default: the genotypes at the 10/50/90 %
#                   quantile of the amplitude a of each shape (3 per panel)
#   N_CURVE_POINTS  integer, default 200   points of the smooth mother curves
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R", "pipeline", "load_pipeline.R"))

############## load

if (!exists("H2")) H2 <- 0.4
if (!exists("CAUSAL_SNP_NUM")) CAUSAL_SNP_NUM <- 5
if (!exists("NUM_GENOTYPES")) NUM_GENOTYPES <- 800
if (!exists("NUM_SNPs")) NUM_SNPs <- 2000
if (!exists("SEED")) SEED <- 1
if (!exists("N_CURVE_POINTS")) N_CURVE_POINTS <- 200

sim <- simulate_reaction_norms(sim_config(h2 = H2, Mc = CAUSAL_SNP_NUM, n_per_shape = NUM_GENOTYPES / 4,
                                          n_snps = NUM_SNPs, seed = SEED))
phenotypes <- sim$phenotypes
mother_par <- sim$mother_parameters
cfg        <- sim$config

shape_label <- c(linear = "Linear", gaussian = "Gaussian", sinusoidal = "Sinusoidal", wave = "Wave")
shape_levels <- c("Gaussian", "Sinusoidal", "Wave", "Linear")

if (!exists("PLOT_GENOTYPES")) {
  PLOT_GENOTYPES <- lapply(setNames(cfg$shapes, cfg$shapes), function(s) {
    mp <- mother_par[mother_par$Shape == s, ]
    mp$Genotype[sapply(quantile(mp$a, c(0.1, 0.5, 0.9)), function(q) which.min(abs(mp$a - q)))]
  })
}

# Measured values (replicates) of all shapes in one table
measured <- bind_rows(lapply(names(phenotypes), function(s) {
  mutate(phenotypes[[s]], ReactionNorm = shape_label[[s]])
})) %>% mutate(ReactionNorm = factor(ReactionNorm, levels = shape_levels))

# Smooth mother reaction norms of the given genotypes, from their mother parameters
mother_curves <- function(shape, ids) {
  mp <- mother_par[mother_par$Shape == shape, ]
  mp <- mp[match(ids, mp$Genotype), ]
  t <- seq(0, 1, length.out = N_CURVE_POINTS)
  Y <- reaction_norm_matrix(shape, t, as.list(mp[reaction_norm_parameters(cfg, shape)]), cfg)
  data.frame(Genotype = rep(ids, each = length(t)), Environment = rep(1 + 9 * t, times = length(ids)),
             Trait = as.vector(t(Y)), ReactionNorm = factor(shape_label[[shape]], levels = shape_levels))
}

############## histograms of the measured trait values

make_hist = function(subdf, title, show_axes = TRUE) {
  ggplot(subdf, aes(x = Trait)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins  = 30,
                   fill  = "grey70",
                   color = "white") +
    geom_density(linewidth = 1, color = "black") +
    labs(x = if (show_axes) "Trait value" else NULL,
         y = if (show_axes) "Density"     else NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.text        = if (show_axes) element_text() else element_blank(),
      axis.ticks       = if (show_axes) element_line() else element_blank()
    )
}

p1 = make_hist(filter(measured, ReactionNorm == "Gaussian"),   "Gaussian")
p2 = make_hist(filter(measured, ReactionNorm == "Sinusoidal"), "Sinusoidal")
p3 = make_hist(filter(measured, ReactionNorm == "Wave"),       "Wave")
p4h = make_hist(filter(measured, ReactionNorm == "Linear"),    "Linear")

grid_plot = (p1 | p2) /
  (p3 | p4h)
print(grid_plot)

p_combined = make_hist(measured, "Combined Genotype Forms", show_axes = TRUE)
print(p_combined)

############## reaction-norm panels

# mother curve (solid), replicate measurements (dashed) and replicate range (ribbon)
# for the genotypes `ids` of every shape; colour by `colour_by` ("ReactionNorm" or "Genotype")
reaction_norm_layers <- function(ids, colour_by, show_ribbon = TRUE, show_dashed = TRUE) {
  shapes <- names(ids)
  mothers <- bind_rows(lapply(shapes, function(s) mother_curves(s, ids[[s]])))
  reps <- filter(measured, Genotype %in% unlist(ids))
  ribbon <- reps %>% group_by(ReactionNorm, Genotype, Environment) %>%
    summarise(ymin = min(Trait), ymax = max(Trait), .groups = "drop")
  list(
    if (show_ribbon) geom_ribbon(data = ribbon, aes(x = Environment, ymin = ymin, ymax = ymax,
                                                    fill = .data[[colour_by]], group = Genotype),
                                 alpha = 0.2, show.legend = FALSE),
    geom_line(data = mothers, aes(x = Environment, y = Trait, color = .data[[colour_by]], group = Genotype),
              linewidth = 1),
    if (show_dashed) geom_line(data = reps, aes(x = Environment, y = Trait, color = .data[[colour_by]],
                                               group = interaction(Genotype, Replicate)),
                               linetype = "dashed", linewidth = 0.5, show.legend = FALSE)
  )
}

line_theme <- theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank()
  )

# one genotype of every shape in one plot
create_plot4 = function(ids, show_ribbon = TRUE, show_dashed = TRUE) {
  ggplot() +
    reaction_norm_layers(lapply(ids, `[`, 1), "ReactionNorm", show_ribbon, show_dashed) +
    labs(title = "", x = "Environment", y = "Trait Value", color = "Reaction-Norm Type", fill = NULL) +
    line_theme +
    theme(
      legend.position      = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      legend.background    = element_rect(fill = scales::alpha("white", 0.5), linewidth = 0),
      legend.key.size      = unit(0.3, "lines"),
      legend.spacing.x     = unit(0.4, "lines"),
      legend.margin        = margin(2, 2, 2, 2, "pt")
    )
}

p4 = create_plot4(PLOT_GENOTYPES)
print(p4)

# one panel per shape with its selected genotypes
create_multi_panel_plot = function(ids, show_ribbon = TRUE, show_dashed = FALSE) {
  plots = lapply(c("gaussian", "sinusoidal", "wave", "linear"), function(s) {
    ggplot() +
      reaction_norm_layers(ids[s], "Genotype", show_ribbon, show_dashed) +
      labs(x = "Environment", y = "Trait Value", title = shape_label[[s]]) +
      line_theme + theme(legend.position = "none", plot.title = element_text(size = 9))
  })
  (plots[[1]] | plots[[2]]) /
    (plots[[3]] | plots[[4]])
}

p = create_multi_panel_plot(PLOT_GENOTYPES, show_ribbon = TRUE, show_dashed = TRUE)
print(p)

################### combined

design <- "
AAAAA
AAAAA
AAAAA
BBBBB
BBBBB
BBBBB
"

final_fig2 <- wrap_plots(
  A = wrap_elements(full = p),           # 2x2 line panels
  B = wrap_elements(full = grid_plot),   # 2x2 histograms
  design = design
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))
print(final_fig2)

final_fig3 <- wrap_plots(
  A = wrap_elements(full = p4),          # one genotype of every shape
  B = wrap_elements(full = grid_plot),   # 2x2 histograms
  design = design
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))
print(final_fig3)

ggsave(create.dir = TRUE,
  filename = here::here("output", "plots", "reaction_norms_final_figure.pdf"),
  plot     = final_fig3,
  device   = "pdf",
  width    = 6.3,
  height   = 7,
  units    = "in",
  dpi      = 900
)
ggsave(create.dir = TRUE,
  filename = here::here("output", "plots", "reaction_norms_final_figure_panels.pdf"),
  plot     = final_fig2,
  device   = "pdf",
  width    = 6.3,
  height   = 7,
  units    = "in",
  dpi      = 900
)
ggsave(create.dir = TRUE,
  filename = here::here("output", "plots", "p_combined.pdf"),
  plot     = p_combined,
  device   = "pdf",
  width    = 6.3,
  height   = 4,
  units    = "in",
  dpi      = 900
)
