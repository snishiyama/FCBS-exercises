library(mvtnorm)

(dat <- scan("data/swim.dat") |>
  matrix(nrow = 4, byrow = TRUE))

(week <- seq(from = 0, by = 2, length.out = 6))

gibbs <- \(y, iv, sig2_0, n_iter = 1000) {
  n <- length(y)
  # prior
  beta_0 <- c(23, 0)
  Sig_0 <- matrix(c(0.5^2, 0, 0, 0.1^2), nrow = 2)
  nu0 <- 2

  # IV
  beta1 <- rep(1, times = length(iv))
  X <- cbind(beta1, iv)

  iSig_0 <- solve(Sig_0)
  XX <- t(X) %*% X

  # init
  Beta <- matrix(nrow = n_iter, ncol = ncol(X))
  sig2 <- NULL

  B <- beta_0
  s2 <- var(y)
  # main loop
  for (i in seq_len(n_iter)) {
    # Beta
    V <- solve(iSig_0 + XX / s2)
    E <- V %*% (iSig_0 %*% beta_0 + t(X) %*% y / s2)
    B <- t(mvtnorm::rmvnorm(1, E, V))

    # sig2
    ssr <- sum(y^2) - 2 * t(B) %*% t(X) %*% y + t(B) %*% t(X) %*% X %*% B
    s2 <- 1 / rgamma(1, (nu0 + n) / 2, (nu0 * sig2_0 + ssr) / 2)

    Beta[i, ] <- B
    sig2 <- c(sig2, s2)
  }

  return(list(beta = Beta, sig2 = sig2))
}


ptys <- c(19, 22, 23, 24) # point types
plot(week, dat[1, ], pch = ptys[1], ylim = c(min(dat), max(dat)))
lines(week, dat[1, ])
for (i in 2:4) {
  lines(week, dat[i, ])
  points(week, dat[i, ], pch = ptys[i])
}

n_iter <- 2000
sig2_0 <- 0.5^2

post_1 <- gibbs(dat[1, ], week, sig2_0, n_iter)
post_2 <- gibbs(dat[2, ], week, sig2_0, n_iter)
post_3 <- gibbs(dat[3, ], week, sig2_0, n_iter)
post_4 <- gibbs(dat[4, ], week, sig2_0, n_iter)

pred <- \(iv, Beta, sig2) {
  const <- rep(1, times = length(iv))
  X <- cbind(const, iv)
  n_iter <- nrow(Beta)
  mat <- matrix(nrow = n_iter, ncol = nrow(X))
  for (i in seq_len(n_iter)) {
    XB <- X %*% Beta[i, ]
    samp <- rmvnorm(1, XB, sig2[i] * diag(nrow(XB)))
    mat[i, ] <- samp
  }
  return(mat)
}


# a-i --------------------------------------------------

pred_1 <- pred(week, post_1$beta, post_1$sig2)
pred_2 <- pred(week, post_2$beta, post_2$sig2)
pred_3 <- pred(week, post_3$beta, post_3$sig2)
pred_4 <- pred(week, post_4$beta, post_4$sig2)


plotter <- \(x, y_dat, y_pred) {
  qts <- apply(y_pred, 2, \(x) quantile(x, c(0.025, 0.5, 0.975)))
  plot(x, y_dat, pch = ptys[1], ylim = c(21.5, 24.5))
  lines(x, qts[1, ], lty = "dashed")
  lines(x, qts[2, ])
  lines(x, qts[3, ], lty = "dashed")
}

par(mfrow = c(2, 2))
preds <- list(pred_1, pred_2, pred_3, pred_4)
for (i in 1:4) {
  plotter(week, dat[i, ], preds[[i]])
}
dev.off()


# a-ii --------------------------------------------------

week_extra <- c(week, 12)
extra_1 <- pred(week_extra, post_1$beta, post_1$sig2)
extra_2 <- pred(week_extra, post_2$beta, post_2$sig2)
extra_3 <- pred(week_extra, post_3$beta, post_3$sig2)
extra_4 <- pred(week_extra, post_4$beta, post_4$sig2)

plotter_extra <- \(x, y_dat, y_pred) {
  qts <- apply(y_pred, 2, \(x) quantile(x, c(0.025, 0.5, 0.975)))
  plot(x[-length(x)], y_dat,
    pch = 19,
    xlim = c(min(x), max(x)), ylim = c(21.5, 24.5),
    xlab = "week", ylab = "time [s]"
  )
  lines(x, qts[1, ], lty = "dashed")
  lines(x, qts[2, ])
  lines(x, qts[3, ], lty = "dashed")
}

par(mfrow = c(2, 2))
extras <- list(extra_1, extra_2, extra_3, extra_4)
for (i in 1:4) {
  plotter_extra(week_extra, dat[i, ], extras[[i]])
}
dev.off()


# b --------------------------------------------------

pred_week12 <- cbind(extra_1[, 7], extra_2[, 7], extra_3[, 7], extra_4[, 7])
table(apply(pred_week12, 1, which.max)) / n_iter
table(apply(pred_week12, 1, which.min)) / n_iter
