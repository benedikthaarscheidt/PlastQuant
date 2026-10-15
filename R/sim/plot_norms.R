# =============================================================================
# plot_norms.R — the true reaction norms of a simulated population
# =============================================================================
# Drawn from the parameters (no measurement): part of simulate_population(), saved
# as data/reaction_norms.pdf next to the norms themselves.

# Smooth true curves of the given genotypes: their mother norm (Replicate 0) and the
# norms of their replicate plants, at n_grid equidistant points of the gradient
true_norm_curves <- function(norms, shape, ids, n_grid = 200) {
  cfg <- norms$config
  thetas <- reaction_norm_parameters(cfg, shape)
  t <- (seq_len(n_grid) - 1) / (n_grid - 1)
  mp <- norms$mother_parameters[norms$mother_parameters$Shape == shape, ]
  mp <- mp[mp$Genotype %in% ids, ]
  pp <- norms$plant_parameters[[shape]]
  pp <- pp[pp$Genotype %in% ids, ]
  Ym <- reaction_norm_matrix(shape, t, as.list(mp[thetas]), cfg)
  Yp <- reaction_norm_matrix(shape, t, as.list(pp[thetas]), cfg)
  rbind(long_table(mp$Genotype, rep(0L, nrow(mp)), 1 + 9 * t, Ym),
        long_table(pp$Genotype, pp$Replicate, 1 + 9 * t, Yp))
}

# One page per shape: for the genotypes at the 10, 50 and 90 % quantile of the
# amplitude a, the mother norm (solid) and the replicate plants (dashed), on a fixed
# trait axis from 0 to 2 with the legend right of the plot. Writes a PDF to file.
plot_true_norms <- function(norms, file, n_grid = 200) {
  cfg <- norms$config
  cols <- c("#1b9e77", "#d95f02", "#7570b3")
  grDevices::pdf(file, width = 8, height = 5)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(5, 4, 4, 11))   # room for the legend right of the plot
  for (shape in cfg$shapes) {
    mp <- norms$mother_parameters[norms$mother_parameters$Shape == shape, ]
    picked <- mp$Genotype[unique(vapply(stats::quantile(mp$a, c(0.1, 0.5, 0.9)),
                                        function(q) which.min(abs(mp$a - q)), integer(1)))]
    sel <- true_norm_curves(norms, shape, picked, n_grid)
    graphics::plot(range(sel$Environment), c(0, 2), type = "n", ylim = c(0, 2),
                   xlab = "Environment", ylab = "Trait value", main = paste0(shape, ": true reaction norms"))
    for (k in seq_along(picked)) {
      for (r in sort(unique(sel$Replicate[sel$Genotype == picked[k]]))) {
        d <- sel[sel$Genotype == picked[k] & sel$Replicate == r, ]
        graphics::lines(d$Environment, d$Trait, col = cols[k], lty = if (r == 0) 1 else 2, lwd = if (r == 0) 2 else 1)
      }
    }
    graphics::legend("topleft", inset = c(1.02, 0), xpd = TRUE, bty = "n", cex = 0.8, col = c(cols[seq_along(picked)], "grey40", "grey40"),
                     lty = c(rep(1, length(picked)), 1, 2), lwd = c(rep(2, length(picked)), 2, 1),
                     legend = c(picked, "mother norm", "replicate plants"))
  }
  invisible(file)
}
