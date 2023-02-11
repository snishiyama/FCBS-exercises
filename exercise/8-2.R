# a --------------------------------------------------

gibbs_sampling <- function(del0, t02, n_iter = 5000) {
  # prior
  mu0 <- 75
  g02 <- 100
  # del0 <-
  # t02 <-
  s20 <- 100
  nu0 <- 2

  # data
  n_a <- n_b <- 16
  ybar_a <- 75.2
  s_a <- 7.3
  ybar_b <- 77.5
  s_b <- 8.1

  sum_a <- ybar_a * n_a
  sum_b <- ybar_b * n_b
  ss_a <- n_a * s_a^2 + n_a * ybar_a^2 # sum of yi_a ^ 2
  ss_b <- n_b * s_b^2 + n_b * ybar_b^2 # sum of yi_b ^ 2

  # gibbs
  # initialize
  mu <- (ybar_a + ybar_b) / 2
  del <- (ybar_a - ybar_b) / 2
  MU <- DEL <- S2 <- NULL
  for (s in seq_len(n_iter)) {
    # s2
    nun <- nu0 + n_a + n_b # nun
    nun_s2n <- nu0 * s20 +
      ss_a - 2 * sum_a * (mu + del) + n_a * (mu + del)^2 +
      ss_b - 2 * sum_b * (mu - del) + n_b * (mu - del)^2
    s2 <- 1 / rgamma(1, nun / 2, nun_s2n / 2)

    # mu
    gn2 <- 1 / (1 / g02 + (n_a + n_b) / s2)
    mun <- gn2 * (mu0 / g02 + (sum_a + sum_b) / s2)
    mu <- rnorm(1, mun, sqrt(gn2))

    # del
    tn2 <- 1 / (1 / t02 + (n_a + n_b) / s2)
    deln <- tn2 * (del0 / t02 + (sum_a - sum_b) / s2)
    del <- rnorm(1, deln, sqrt(tn2))

    MU <- c(MU, mu)
    DEL <- c(DEL, del)
    S2 <- c(S2, s2)
  }

  return(list(MU = MU, DEL = DEL, S2 = S2))
}

opts_del0 <- seq(-4, 4, by = 2)
opts_t02 <- c(10, 50, 100, 500)
seq_del0 <- rep(opts_del0, times = length(opts_t02))
seq_t02 <- rep(opts_t02, each = length(opts_del0))
(comb_del0_t02 <- rbind(seq_del0, seq_t02)) # just for check

post_samples <- mapply(gibbs_sampling, del0 = seq_del0, t02 = seq_t02)
str(post_samples)


# i --------------------------------------------------
prob_under0 <- sapply(post_samples["DEL", ], \(x) mean(x < 0))
filled.contour(opts_del0, opts_t02, matrix(prob_under0, nrow = 5))


# ii --------------------------------------------------

cis <- lapply(post_samples["DEL", ], \(x) quantile(x, c(0.025, 0.975)))

# png("fig/8-2_a-ii.png", width = 1500, height = 1200, res = 100)
par(mfrow = c(4, 5), mar = c(3, 3, 4, 1))
for (i in seq_along(seq_del0)) {
  title <- bquote(delta[0] == .(seq_del0[i]) * "," ~ tau[0]^2 == .(seq_t02[i]))
  plot(density(post_samples["DEL", i][[1]], adjust = 2),
    xlim = c(-6, 6), ylim = c(0, .32),
    lwd = 2, main = title
  )
  abline(v = cis[[i]][1]) # lower
  abline(v = cis[[i]][2]) # upper
  txt <- bquote(Pr(delta < 0 ~ "|" ~ Y) == .(sprintf("%.2f", prob_under0[i])))
  legend("topright", legend = txt, x.intersp = -0.5)
}
dev.off()


# iii --------------------------------------------------

mc_prior <- \(del0, t02, n = 5000) {
  MU <- rnorm(n, 75, 100)
  DEL <- rnorm(n, del0, sqrt(t02))
  return(list(MU = MU, DEL = DEL))
}

calc_corr <- \(mu, del) {
  cor(mu + del, mu - del)
}

pri_samples <- mapply(mc_prior, del0 = seq_del0, t02 = seq_t02)
pri_corrs <- apply(pri_samples, 2, \(x) calc_corr(x["MU"][[1]], x["DEL"][[1]]))

post_corrs <- apply(
  post_samples, 2,
  \(x) calc_corr(x["MU"][[1]], x["DEL"][[1]])
)

filled.contour(opts_del0, opts_t02, matrix(pri_corrs, nrow = 5))
filled.contour(opts_del0, opts_t02, matrix(post_corrs, nrow = 5))

# png("fig/8-2_a-iii.png", width = 1200, height = 600, res = 100)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 1))
contour(opts_del0, opts_t02, matrix(pri_corrs, nrow = 5),
  labcex = 1.2,
  xlab = expression(delta[0]), ylab = expression(tau[0])
)
contour(opts_del0, opts_t02, matrix(post_corrs, nrow = 5),
  labcex = 1.2,
  xlab = expression(delta[0]), ylab = expression(tau[0])
)
dev.off()
