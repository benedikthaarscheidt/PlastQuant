# =============================================================================
# plot_captions.R — captions under the figures of the analysis scripts
# =============================================================================
# Every figure gets a short caption in two parts: what is shown (the setup) and what
# can be read from it. The captions are added only where a figure is printed or saved,
# so composite figures built from the same panels get one caption, not one per panel.
# options(plastquant.captions = FALSE) switches them off (e.g. for the paper figures,
# whose legends are in the manuscript).

captions_on <- function() isTRUE(getOption("plastquant.captions", TRUE))

# Caption text: "<shown> <reading>", wrapped to `width` characters
caption_text <- function(shown, reading, width = 120) {
  paste(strwrap(paste(shown, reading), width = width), collapse = "\n")
}

caption_theme <- function() {
  ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7, colour = "grey25"))
}

# ggplot or patchwork: the caption under the (whole) figure; existing annotations such
# as panel tags are kept
with_caption <- function(p, shown, reading, width = 120) {
  if (!captions_on()) return(p)
  txt <- caption_text(shown, reading, width)
  if (inherits(p, "patchwork")) {
    p + patchwork::plot_annotation(caption = txt, theme = caption_theme())
  } else {
    p + ggplot2::labs(caption = txt) + caption_theme()
  }
}

# grid objects (gtable, gridExtra arrangements, pheatmap$gtable): caption as a text row below
grob_with_caption <- function(g, shown, reading, width = 120) {
  if (!captions_on()) return(g)
  gridExtra::arrangeGrob(g, bottom = caption_grob(shown, reading, width))
}

# Caption as a grob, e.g. for grid.arrange(..., bottom = caption_grob(...)); NULL when
# captions are off (grid.arrange then draws no caption row)
caption_grob <- function(shown, reading, width = 120) {
  if (!captions_on()) return(NULL)
  grid::textGrob(caption_text(shown, reading, width), x = grid::unit(0.01, "npc"), hjust = 0,
                 gp = grid::gpar(fontsize = 7, col = "grey25"))
}

# Base graphics: caption in the outer bottom margin of the current page. Call
# caption_margin() before the plot (it reserves the margin) and base_caption() after it.
caption_margin <- function(lines = 4) {
  if (captions_on()) graphics::par(oma = c(lines, 0, 0, 0))
}
base_caption <- function(shown, reading, width = 120) {
  if (!captions_on()) return(invisible(NULL))
  graphics::mtext(caption_text(shown, reading, width), side = 1, outer = TRUE, adj = 0,
                  line = 0.5, padj = 0, cex = 0.6, col = "grey25")
}

# Base graphics with par(mfrow): caption once per page. Call after every panel with
# end = FALSE (draws when the page is full) and once after the last panel with end = TRUE
# (draws when that page is not full, i.e. has no caption yet).
base_caption_panel <- function(shown, reading, end = FALSE, width = 120) {
  mfg <- graphics::par("mfg")
  full <- mfg[1] == mfg[3] && mfg[2] == mfg[4]
  if (full != end) base_caption(shown, reading, width)
}

# Grid page (e.g. a table drawn with grid.draw): caption at the bottom of the page
grid_page_caption <- function(shown, reading, width = 120) {
  g <- caption_grob(shown, reading, width)
  if (is.null(g)) return(invisible(NULL))
  g$y <- grid::unit(0.01, "npc"); g$vjust <- 0
  grid::grid.draw(g)
}

# Setup sentence of a resolution study: population and the resolutions compared
resolutions_setup <- function(population_dir, n_points) {
  paste0("Population ", basename(normalizePath(population_dir)), "; the same plants measured at n = ",
         paste(n_points, collapse = ", "), " equidistant points; PPIs of the replicate-averaged reaction norm of each genotype.")
}

# Setup sentence of a resolution folder: population, measurement points, replicate
# handling of the PPIs and, if given, the genotypes per shape (named counts)
resolution_setup <- function(resolution_dir, handling = NULL, genotypes = NULL) {
  pop <- basename(dirname(normalizePath(resolution_dir)))
  n <- sub("^n_", "", basename(resolution_dir))
  paste0(
    "Population ", pop, ", reaction norms measured at ", n, " points",
    if (!is.null(handling)) paste0(", PPIs computed from ",
      if (handling == "mean_before") "the replicate-averaged reaction norm of each genotype" else "each replicate plant"),
    if (!is.null(genotypes)) paste0("; genotypes: ", paste(genotypes, names(genotypes), collapse = ", ")),
    "."
  )
}
