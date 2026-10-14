# =============================================================================
# plotting_reaction_norms_all.R — reaction-norm overview figure of a population
# =============================================================================
# WHAT IT DOES: Draws the reaction norms of a population folder: per shape (or trait)
#   selected genotypes with their replicate measurements (dashed) and replicate range
#   (ribbon), plus histograms of the measured trait values. For a simulated population
#   the solid line is the mother reaction norm (smooth curve from the mother
#   parameters); for own data it is the genotype mean over the replicates.
# REQUIRES:     A population folder (data/), simulated with simulate_population() /
#               a scenario driver, or own data placed in data/.
# PRODUCES:     POPULATION_DIR/plots/reaction_norms_final_figure.pdf,
#               reaction_norms_final_figure_panels.pdf, trait_histogram_combined.pdf.
# HOW TO RUN:   setwd("~/PlastQuant")
#               POPULATION_DIR <- here::here("output", "populations", "my_population")
#               source(here::here("R","figures","plotting_reaction_norms_all.R"))
# -----------------------------------------------------------------------------
# PARAMETERS — set these before sourcing (guarded defaults)
#   POPULATION_DIR  path, REQUIRED   population folder to draw                            [COMMON]
#   PLOT_GENOTYPES  list shape/trait -> genotype IDs, default: 3 genotypes per panel at the
#                   10/50/90 % quantile of the amplitude a (simulated) or of the genotype
#                   mean trait value (own data)
#   N_CURVE_POINTS  integer, default 200   points of the smooth mother curves (simulated)
# =============================================================================
options(warn = -1)  # silence warnings even if the project .Rprofile was not loaded

library(ggplot2)
library(dplyr)
library(patchwork)
source(here::here("R", "pipeline", "load_pipeline.R"))

############## load

if (!exists("POPULATION_DIR")) stop("set POPULATION_DIR (a population folder) first")
if (!exists("N_CURVE_POINTS")) N_CURVE_POINTS <- 200

pop        <- read_population(POPULATION_DIR)
phenotypes <- pop$phenotypes
simulated  <- !is.null(pop$config)
cfg        <- pop$config
mother_par <- pop$mother_parameters

# Panel order and labels: the shapes of the model in the order of the manuscript figure,
# or the traits of own data as they are
shape_label <- c(linear = "Linear", gaussian = "Gaussian", sinusoidal = "Sinusoidal", wave = "Wave")
panels <- if (simulated) intersect(c("gaussian", "sinusoidal", "wave", "linear"), names(phenotypes)) else names(phenotypes)
labels <- setNames(ifelse(panels %in% names(shape_label), shape_label[panels], panels), panels)

if (!exists("PLOT_GENOTYPES")) {
  pick3 <- function(ids, value) ids[unique(sapply(quantile(value, c(0.1, 0.5, 0.9)), function(q) which.min(abs(value - q))))]
  PLOT_GENOTYPES <- lapply(setNames(panels, panels), function(s) {
    if (simulated) {
      mp <- mother_par[mother_par$Shape == s, ]
      pick3(mp$Genotype, mp$a)
    } else {
      g <- tapply(phenotypes[[s]]$Trait, phenotypes[[s]]$Genotype, mean, na.rm = TRUE)
      pick3(names(g), as.numeric(g))
    }
  })
}

# Measured values (replicates) of all panels in one table
measured <- bind_rows(lapply(panels, function(s) {
  mutate(phenotypes[[s]], ReactionNorm = labels[[s]])
})) %>% mutate(ReactionNorm = factor(ReactionNorm, levels = unname(labels)))

# Solid line of the given genotypes: smooth mother reaction norm (simulated) or the
# genotype mean over the replicates (own data)
center_curves <- function(shape, ids) {
  if (simulated) {
    mp <- mother_par[mother_par$Shape == shape, ]
    mp <- mp[match(ids, mp$Genotype), ]
    t <- seq(0, 1, length.out = N_CURVE_POINTS)
    Y <- reaction_norm_matrix(shape, t, as.list(mp[reaction_norm_parameters(cfg, shape)]), cfg)
    out <- data.frame(Genotype = rep(ids, each = length(t)), Environment = rep(1 + 9 * t, times = length(ids)),
                      Trait = as.vector(t(Y)))
  } else {
    out <- phenotypes[[shape]] %>% filter(Genotype %in% ids) %>%
      group_by(Genotype, Environment) %>% summarise(Trait = mean(Trait, na.rm = TRUE), .groups = "drop")
  }
  mutate(as.data.frame(out), ReactionNorm = factor(labels[[shape]], levels = unname(labels)))
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
         y = if (show_axes) "Density"     else NULL,
         title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      axis.text        = if (show_axes) element_text() else element_blank(),
      axis.ticks       = if (show_axes) element_line() else element_blank()
    )
}

grid_plot = wrap_plots(lapply(panels, function(s) make_hist(filter(measured, ReactionNorm == labels[[s]]), labels[[s]])),
                       ncol = 2)
print(grid_plot)

p_combined = make_hist(measured, "All shapes/traits", show_axes = TRUE)
print(p_combined)

############## reaction-norm panels

# solid line, replicate measurements (dashed) and replicate range (ribbon) for the
# genotypes `ids` of every panel; colour by `colour_by` ("ReactionNorm" or "Genotype")
reaction_norm_layers <- function(ids, colour_by, show_ribbon = TRUE, show_dashed = TRUE) {
  centers <- bind_rows(lapply(names(ids), function(s) center_curves(s, ids[[s]])))
  reps <- bind_rows(lapply(names(ids), function(s) {
    filter(measured, ReactionNorm == labels[[s]], Genotype %in% ids[[s]])
  }))
  ribbon <- reps %>% group_by(ReactionNorm, Genotype, Environment) %>%
    summarise(ymin = min(Trait), ymax = max(Trait), .groups = "drop")
  list(
    if (show_ribbon) geom_ribbon(data = ribbon, aes(x = Environment, ymin = ymin, ymax = ymax,
                                                    fill = .data[[colour_by]], group = Genotype),
                                 alpha = 0.2, show.legend = FALSE),
    geom_line(data = centers, aes(x = Environment, y = Trait, color = .data[[colour_by]], group = Genotype),
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
    axis.ticks  = element_blank(),
    panel.grid  = element_blank()
  )

# one genotype of every panel in one plot
create_plot4 = function(ids, show_ribbon = TRUE, show_dashed = TRUE) {
  ggplot() +
    reaction_norm_layers(lapply(ids, `[`, 1), "ReactionNorm", show_ribbon, show_dashed) +
    labs(title = "", x = "Environment", y = "Trait Value", color = NULL, fill = NULL) +
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

# one panel per shape/trait with its selected genotypes
create_multi_panel_plot = function(ids, show_ribbon = TRUE, show_dashed = FALSE) {
  plots = lapply(panels, function(s) {
    ggplot() +
      reaction_norm_layers(ids[s], "Genotype", show_ribbon, show_dashed) +
      labs(x = "Environment", y = "Trait Value", title = labels[[s]]) +
      line_theme + theme(legend.position = "none", plot.title = element_text(size = 9))
  })
  wrap_plots(plots, ncol = 2)
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
  A = wrap_elements(full = p),           # line panels
  B = wrap_elements(full = grid_plot),   # histograms
  design = design
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))
print(final_fig2)

final_fig3 <- wrap_plots(
  A = wrap_elements(full = p4),          # one genotype of every panel
  B = wrap_elements(full = grid_plot),   # histograms
  design = design
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))
print(final_fig3)

plots_dir <- population_subdir(POPULATION_DIR, "plots")
ggsave(filename = file.path(plots_dir, "reaction_norms_final_figure.pdf"), plot = final_fig3,
       device = "pdf", width = 6.3, height = 7, units = "in", dpi = 900)
ggsave(filename = file.path(plots_dir, "reaction_norms_final_figure_panels.pdf"), plot = final_fig2,
       device = "pdf", width = 6.3, height = 7, units = "in", dpi = 900)
ggsave(filename = file.path(plots_dir, "trait_histogram_combined.pdf"), plot = p_combined,
       device = "pdf", width = 6.3, height = 4, units = "in", dpi = 900)
message("Plots written to ", plots_dir)
