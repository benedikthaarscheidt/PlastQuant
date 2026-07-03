#this is the script for the assembly of the first figure containing the genotype line plots as well as the histograms for the trait distribution 


library(readr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dplyr)


############## load and subset 

all_combined_data <- read_csv(here::here("output", "scenario_20_5", "synthetic_data", "fixed_full", "all_combined_data.csv"))
gaussian=all_combined_data[all_combined_data[,5]=="Gaussian",]
wave=all_combined_data[all_combined_data[,5]=="Wave",]
sinusoid=all_combined_data[all_combined_data[,5]=="Sinusoidal",]


wave=as.data.frame(wave)
xvals <- as.numeric(wave[1001:1050, "Environment"])
yvals <- as.numeric(wave[1151:1200, 4])
xlim <- range(xvals)
yr   <- range(yvals)
ylim <- c( min(0, yr[1]) , yr[2] )  

plot(
  x    = xvals,
  y    = yvals,
  type = "l",
  xlab = "Trait 3",
  ylab = "Trait 4",
  xlim = xlim,
  ylim = ylim
)
###############
# 1) Load & tag

linear <- read_csv(here::here("output", "scenario_20_5", "synthetic_data", "fixed_full", "linear_reaction_norms_data.csv")) %>%
  mutate(ReactionNorm = "Linear")

df = bind_rows(all_combined_data, linear) %>%
  mutate(ReactionNorm = factor(ReactionNorm,
                               levels = c("Gaussian", "Sinusoidal", "Wave", "Linear")
  ))

# 2) Helper for a single grey‐fill histogram + black density
make_hist = function(subdf, title, show_axes = TRUE) {
  ggplot(subdf, aes(x = Trait)) +
    geom_histogram(aes(y = ..density..),
                   bins  = 30,
                   fill  = "grey70",
                   color = "white") +
    geom_density(size = 1, color = "black") +
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

# 3) Build the four separate plots (no legends needed)
p1 = make_hist(filter(df, ReactionNorm == "Gaussian"),   "Gaussian")
p2 = make_hist(filter(df, ReactionNorm == "Sinusoidal"), "Sinusoidal")
p3 = make_hist(filter(df, ReactionNorm == "Wave"),       "Wave")
p4 = make_hist(filter(df, ReactionNorm == "Linear"),     "Linear")

# 4) Display as 2×2
grid_plot = (p1 | p2) /
  (p3 | p4)
print(grid_plot)

# 5) Combined histogram of all data
p_combined = make_hist(df, "Combined Genotype Forms", show_axes = TRUE)
print(p_combined)

############## plotting 



get_genotype_at_row = function(df, n) {
  unique(df$Genotype)[
    ceiling(n / ((1 + 3) * length(unique(df$Environment))))
  ]
}


n=1 # this is for the selection of the block of genotypes to plot

create_plot4 = function(gaussian_data, sinusoidal_data, wave_data, linear_data,
                        n,
                        show_ribbon = TRUE,
                        show_dashed = T) {
  # — helper to find the genotype ID block at row n —
  get_genotype_at_row = function(df, n) {
    # each genotype block = 1 mother + 3 replicates × length(env)
    block_size = (1 + 3) * length(unique(df$Environment))
    idx = ceiling(n / block_size)
    unique(df$Genotype)[idx]
  }
  
  # pick one genotype from each dataset
  id_gauss = get_genotype_at_row(gaussian_data, n+200)
  id_sine  = get_genotype_at_row(sinusoidal_data, n+800)
  id_wave  = get_genotype_at_row(wave_data,      n+1200)
  id_lin   = get_genotype_at_row(linear_data,    n+2800)
  
  # bind the four blocks together
  plot_df = bind_rows(
    filter(gaussian_data,   Genotype == id_gauss),
    filter(sinusoidal_data, Genotype == id_sine),
    filter(wave_data,       Genotype == id_wave),
    filter(linear_data,     Genotype == id_lin)
  )
  
  fixed_df     = filter(plot_df, Replicate == 0)
  replicate_df = filter(plot_df, Replicate != 0)
  
  # prepare ribbon data if requested
  if (show_ribbon) {
    ribbon_df = replicate_df %>%
      group_by(ReactionNorm, Environment) %>%
      summarise(
        ymin = min(Trait),
        ymax = max(Trait),
        .groups = "drop"
      )
  }
  
  p = ggplot() +
    # 1) optional ribbon
    { if (show_ribbon) geom_ribbon(
      data = ribbon_df,
      aes(x = Environment, ymin = ymin, ymax = ymax, fill = ReactionNorm),
      alpha = 0.2,
      show.legend = FALSE
    ) } +
    # 2) mother curves
    geom_line(
      data = fixed_df,
      aes(x = Environment, y = Trait, color = ReactionNorm, group = ReactionNorm),
      size = 1
    ) +
    # 3) optional dashed replicates
    { if (show_dashed) geom_line(
      data = replicate_df,
      aes(x = Environment, y = Trait,
          group = interaction(ReactionNorm, Genotype, Replicate),
          color = ReactionNorm),
      linetype = "dashed",
      size = 0.5,
      show.legend = FALSE
    ) } +
    labs(
      title = "",
      x     = "Environment",
      y     = "Trait Value",
      color = "Reaction-Norm Type",
      fill  = NULL
    ) +
    theme_minimal() +
    theme(
      axis.line.x          = element_line(color = "black"),
      axis.line.y          = element_line(color = "black"),
      axis.text.x          = element_blank(),
      axis.ticks           = element_blank(),
      panel.grid           = element_blank(),
      legend.position      = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      legend.background    = element_rect(fill = scales::alpha("white", 0.5), size = 0),
      legend.key.size      = unit(0.3, "lines"),
      legend.spacing.x     = unit(0.4, "lines"),
      legend.margin        = margin(2, 2, 2, 2, "pt")
    )
  
  return(p)
}



p4 = create_plot4(
  gaussian_data    = gaussian,
  sinusoidal_data  = sinusoid,
  wave_data        = wave,
  linear_data      = linear,
  n                = n
)
# display
print(p4)




create_multi_panel_plot = function(gaussian_data,
                                   sinusoidal_data,
                                   wave_data,
                                   linear_data,
                                   idxs,
                                   show_ribbon = TRUE,
                                   show_dashed = FALSE) {
  # helper: from row-index n to actual Genotype ID
  get_genotype_at_row = function(df, n) {
    block_size = (1 + 3) * 50
    b = ceiling(n / block_size)
    unique(df$Genotype)[b]
  }
  
  # build a list of the four data & name pairs
  all_data = list(
    Gaussian   = gaussian_data,
    Sinusoidal = sinusoidal_data,
    Wave       = wave_data,
    Linear     = linear_data
  )
  
  # for each type, pull out the 3 chosen genotypes
  plots = lapply(names(all_data), function(rn) {
    df  = all_data[[rn]]
    ids = sapply(idxs[[rn]], function(n) get_genotype_at_row(df, n))
    sel = df %>% filter(Genotype %in% ids)
    
    fixed_df     = sel %>% filter(Replicate == 0)
    replicate_df = sel %>% filter(Replicate != 0)
    
    if (show_ribbon) {
      ribbon_df = replicate_df %>%
        group_by(Genotype, Environment) %>%
        summarise(
          ymin = min(Trait),
          ymax = max(Trait),
          .groups = "drop"
        )
    }
    
    p = ggplot()
    
    if (show_ribbon) {
      p = p + geom_ribbon(
        data = ribbon_df,
        aes(x = Environment, ymin = ymin, ymax = ymax, fill = factor(Genotype)),
        alpha = 0.2
      )
    }
    
    if (show_dashed) {
      p = p + geom_line(
        data = replicate_df,
        aes(x = Environment, y = Trait,
            group = interaction(Genotype, Replicate),
            color = factor(Genotype)),
        linetype = "dashed", size = 0.4,
        show.legend = FALSE
      )
    }
    
    p = p + geom_line(
      data = fixed_df,
      aes(x = Environment, y = Trait,
          color = factor(Genotype), group = Genotype),
      size = 1
    ) +
      labs(
        x = "Environment",
        y = "Trait Value"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.line.x     = element_line(color = "black"),
        axis.line.y     = element_line(color = "black"),
        axis.text.x     = element_blank(),
        axis.ticks      = element_blank(),
        panel.grid      = element_blank(),
        legend.position = "none"
      )
    
    return(p)
  })
  
  # now stitch them into a 2×2
  combined = (plots[[1]] | plots[[2]]) /
    (plots[[3]] | plots[[4]]) +
    plot_layout(guides = "collect")
  
  return(combined)
}

idxs = list(
  Gaussian   = c(1, 1001, 3201),
  Sinusoidal = c(1, 801, 2401),
  Wave       = c(201, 2801, 1201),
  Linear     = c(1, 201, 401)
)

p = create_multi_panel_plot(
  gaussian_data   = gaussian,
  sinusoidal_data = sinusoid,
  wave_data       = wave,
  linear_data     = linear,
  idxs            = idxs,
  show_ribbon     = TRUE,
  show_dashed     = T
)

print(p)


################### combined



final_fig <- (p) /
  (grid_plot) /
  (p_combined) +
  plot_layout(
    ncol    = 2, 
    heights = c(1, 1.5, 1)   # give the line‐plots more vertical room
  ) +
  plot_annotation()

print(final_fig)

altA <- (p / grid_plot) | p_combined +
  plot_layout(
    heights = c(2, 1),  # left‐col: line panels twice as tall as hist grid
    widths  = c(2, 1)   # left‐col twice as wide as right‐col
  ) +
  plot_annotation()

print(altA)

altB <-(p | grid_plot) /
  p_combined +
  plot_layout(
    heights = c(2, 1)  # give the panel‐row twice the vertical space
  ) +
  plot_annotation()

print(altB)

design <- "
AAAAA
AAAAA
AAAAA
BBBBB
BBBBB
BBBBB

"

# 2) Wrap each composite in wrap_elements()
final_fig2 <- wrap_plots(
  A = wrap_elements(full = p),           # your 2×2 line‐panel grid
  B = wrap_elements(full = grid_plot),   # your 2×2 hist‐grid
  #C = wrap_elements(full = p_combined),  # your combined histogram
  design = design
) +
  # 3) Auto‐tag them A, B, C and add a title
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))

print(final_fig2)



design2 <- "
AAAAA
AAAAA
AAAAA
BBBBB
BBBBB
BBBBB

"

# 2) Wrap your three composites—now A = p4
final_fig3 <- wrap_plots(
  A = wrap_elements(full = p4),           # the “combined” line‐plot across 4 types
  B = wrap_elements(full = grid_plot),    # your 2×2 histograms
  #C = wrap_elements(full = p_combined),   # the single combined histogram
  design = design2
) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.tag = element_text(size = 8, face = "bold"))
  ) &
  theme(plot.tag.position = c(0, 1))

print(final_fig3)

ggsave(
  filename = here::here("output", "plots", "figures", "reaction_norms_final_figure.pdf"),
  plot     = final_fig3,
  device   = "pdf",
  width    = 6.3,
  height   = 7,
  units    = "in",
  dpi      = 900
)
ggsave(
  filename = here::here("output", "plots", "figures", "reaction_norms_final_figure_panels.pdf"),
  plot     = final_fig2,
  device   = "pdf",
  width    = 6.3,
  height   = 7,
  units    = "in",
  dpi      = 900
)

ggsave(
  filename = here::here("output", "plots", "figures", "p_combined.pdf"),
  plot     = p_combined,
  device   = "pdf",
  width    = 6.3,
  height   = 4,
  units    = "in",
  dpi      = 900
)