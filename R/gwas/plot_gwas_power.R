# =============================================================================
# plot_gwas_power.R — heatmaps of the GWAS power per shape
# =============================================================================
# Part of run_population_gwas(): drawn whenever the population has a truth table.

# power: output of evaluate_power() (Score, Parameter, True_SNPs, Found, Power, Shape).
# One overview page with all shapes, then one page per shape with its false positives.
plot_gwas_power <- function(power, file, subtitle = "") {
  heatmap <- function(df, title, sub) {
    df <- df[df$Parameter != "False_Positives", ]
    ggplot2::ggplot(df, ggplot2::aes(x = Parameter, y = Score, fill = Power)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), na.value = "grey90") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(Power), "", sprintf("%.2f", Power))), color = "black", size = 3) +
      ggplot2::labs(title = title, subtitle = sub, x = "Causal parameter (truth of this shape)", y = "Plasticity index") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 10, height = 8)
  on.exit(grDevices::dev.off())
  print(heatmap(power, "GWAS power per shape", subtitle) +
          ggplot2::facet_grid(. ~ Shape, scales = "free_x", space = "free_x"))
  for (s in unique(power$Shape)) {
    fp <- power[power$Shape == s & power$Parameter == "False_Positives", ]
    print(heatmap(power[power$Shape == s, ], paste0("GWAS power: ", s),
                  paste0(subtitle, "\nFalse positives (all PPIs): ", sum(fp$Found))))
  }
  invisible(file)
}
