test_that("latent SNP share equals h2 and measurement error is 5 %", {
  sim <- simulate_reaction_norms(sim_config(h2 = 0.4, Mc = 10, n_per_shape = 2000, n_rep = 2, n_snps = 2000, seed = 5))
  vr <- variance_report(sim, n_mc = 500)
  expect_equal(mean(vr$latent$snp_share), 0.4, tolerance = 0.02)
  expect_true(all(abs(vr$latent$snp_share - 0.4) < 0.12))
  expect_equal(vr$trait$Share[vr$trait$Component == "measurement error"], rep(0.05, 4), tolerance = 0.01)
  expect_equal(as.vector(tapply(vr$trait$Share, vr$trait$Shape, sum)), rep(1, 4), tolerance = 1e-12)
  expect_true(all(vr$measurement$n_negative == 0))
})

test_that("for the linear shape the trait-level SNP share is close to h2", {
  set.seed(9)
  d <- trait_variance_decomposition(sim_config(h2 = 0.6, Mc = 5), "linear", n_mc = 4000)
  expect_equal(d$snp / (d$snp + d$background + d$interaction), 0.6, tolerance = 0.05)
})
