library(coda)
library(mvtnorm)

dat <- read.table("data/prayer.dat", header = TRUE)
str(dat)

# functions for gibbs sampling ------------------------------------------------

# sample from truncated normal
rtnorm <- \(n, mu, sigma, a, b) {
  bound <- \(x) pnorm(x, mu, sigma) # equal to pnorm((x- mu) / sigma, 0, 1)
  u <- runif(n, bound(a), bound(b))
  return(qnorm(u, mu, sigma))
}

sample_beta <- \(z, X) {
  n <- nrow(X)
  vb <- n / (n + 1) * solve(t(X) %*% X)
  eb <- vb %*% (t(X) %*% z)
  return(t(chol(vb)) %*% rnorm(ncol(X)) + eb) # faster than mvtnorm::rmvnorm(1, eb, vb)
}

sample_z_by_rank <- \(z, beta, y, X) {
  ez <- X %*% beta
  zp <- numeric(length(y))
  # y is already rank variables
  for (r in sample(unique(y))) {
    ir <- seq_along(y)[y == r]
    is_smaller <- y < r
    is_larger <- r < y
    a <- ifelse(sum(is_smaller) > 0, max(z[is_smaller]), -Inf)
    b <- ifelse(sum(is_larger) > 0, min(z[is_larger]), Inf)
    zp[ir] <- rtnorm(length(ir), ez[ir], 1, a, b)
  }
  return(zp)
}

gibbs_ranklike <- \(y, X, n_iter) {
  params <- list(
    BETA = matrix(nrow = n_iter, ncol = ncol(X)),
    Z = matrix(nrow = n_iter, ncol = length(y)),
    ACS = 0
  )

  # not sure if this runs faster than sample_beta(z, X)
  sample_beta_ <- \(z) sample_beta(z, X)

  # init
  beta <- rep(0, ncol(X))
  z <- qnorm(rank(y, ties.method = "random") / (length(y) + 1))

  start <- Sys.time()
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) {
      elapsed <- sprintf("%.2f", difftime(Sys.time(), start, unit = "secs"))
      cat("Running", i, "th loop. Elapsed", elapsed, "sec\n")
    }
    z <- sample_z_by_rank(z, beta, y, X)
    beta <- sample_beta_(z)

    # help mixing
    ez <- X %*% beta
    zp <- z + rnorm(1, 0, length(y)^(-1 / 3))
    lhr <- sum(dnorm(zp, ez, 1, log = TRUE) - dnorm(z, ez, 1, log = TRUE))
    if (log(runif(1)) < lhr) {
      z <- zp
      params$ACS <- params$ACS + 1
    }

    params$BETA[i, ] <- beta
    params$Z[i, ] <- z
  }
  return(params)
}

# run --------------------------------------------------

X <- as.matrix(dat[, c("female", "vocab", "female.vocab")])
y <- dat$prayer

n_iter <- 25000
res <- gibbs_ranklike(y, X, n_iter)
res$ACS / n_iter

res_thin <- (\() {
  lag <- 25
  thinning <- seq(lag, n_iter, by = lag)
  return(list(
    BETA = res$BETA[thinning, ],
    Z = res$Z[thinning, ]
  ))
})()

par(mfrow = c(3, 1))
for (i in 1:3) {
  plot(res_thin$BETA[, i], type = "l")
}
dev.off()

apply(res_thin$BETA, 2, coda::effectiveSize)
summary(apply(res_thin$Z, 2, coda::effectiveSize))

nXtX <- nrow(X) * solve(t(X) %*% X)
BETA_prior <- mvtnorm::rmvnorm(1000, sigma = nXtX)
png("fig/12-1_beta.png", width = 800, height = 400, res = 100)
par(mfrow = c(1, 3), mar = c(3, 3, 3, 1))
for (i in 1:3) {
  plot(density(res_thin$BETA[, i], adj = 2), main = bquote(beta[.(i)]))
  lines(density(BETA_prior[, i], adj = 2), lty = "dashed")
  legend("topright",
    legend = c("Prior", "Posterior"),
    lty = c("dashed", "solid")
  )
}
dev.off()

# scatter plot
eB <- apply(res$BETA, 2, mean)
is_female <- X[, 1] == 1
png("fig/12-1_scatter.png", width = 800, height = 600, res = 100)
plot(X[is_female, 2] + 0.1, res$Z[1000, is_female],
  ylim = range(res$Z[1000, ]),
  xlab = "vocab", ylab = "Z[1000, ]", pch = 21, col = "blue"
)
points(X[!is_female, 2] - 0.1, res$Z[1000, !is_female], pch = 24)
abline(0, eB[2], lwd = 2)
abline(eB[1], eB[2] + eB[3], lwd = 2, col = "blue")
legend(
  x = mean(range(X[, 2])), y = 3.2, xjust = 0.5, yjust = 0,
  legend = c("Male", "Female"), lty = "solid", col = c("black", "blue"),
  xpd = TRUE, horiz = TRUE
)
dev.off()
