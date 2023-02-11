dat <- read.table("data/azdiabetes.dat", header = TRUE)
head(dat)

y <- as.numeric(dat$diabetes == "Yes")
X <- as.matrix(dat[, c("npreg", "bp", "bmi", "ped", "age")])
X <- apply(X, 2, \(x) {
  n <- length(x)
  (x - mean(x)) / sqrt(var(x) * (n - 1) / n)
})
X <- cbind(rep(1, nrow(X)), X)
colnames(X)[1] <- "intercept"
head(X)

# ML
summary(glm(y ~ -1 + X, family = binomial))

# no model selection --------------------------------------------------

metropolis <- \(y, X, sigp, n_iter) {
  p <- ncol(X)
  beta <- rep(0, p)
  acs <- 0
  BETA <- matrix(nrow = n_iter, ncol = p)

  log_prior <- \(beta_) dnorm(beta_, sd = c(16, rep(4, p - 1)))
  log_like <- \(beta_) y * (X %*% beta_) - log(1 + exp(X %*% beta_))

  for (s in seq_len(n_iter)) {
    beta_p <- t(chol(sigp^2 * diag(p))) %*% rnorm(p) + beta
    # beta_p <- t(mvtnorm::rmvnorm(1, beta, sigp^2 * diag(p)))
    lhr <- sum(log_like(beta_p)) - sum(log_like(beta)) +
      sum(log_prior(beta_p)) - sum(log_prior(beta))

    if (log(runif(1)) < lhr) {
      beta <- beta_p
      acs <- acs + 1
    }

    BETA[s, ] <- beta
  }
  return(list(acs = acs, BETA = BETA))
}

n_iter <- 25000
res <- metropolis(y, X, 0.07, n_iter)
res$acs / n_iter
thin <- seq(25, n_iter, by = 25)
acf(res$BETA[thin, ])

par(mfrow = c(2, 3))
for (c in seq_len(ncol(X))) {
  plot(res$BETA[thin, c], type = "l")
}
dev.off()
apply(res$BETA, 2, quantile, c(0.025, 0.5, 0.975))


# a --------------------------------------------------

metro_hasting <- \(y, X, sigp, n_iter) {
  p <- ncol(X)
  beta <- rep(0, p)
  z <- rep(1, p)
  acs <- 0
  BETA <- Z <- matrix(nrow = n_iter, ncol = p)
  colnames(BETA) <- colnames(Z) <- colnames(X)

  log_prior <- \(beta_) dnorm(beta_, sd = c(16, rep(4, p - 1)))
  log_like <- \(beta_) y * (X %*% beta_) - log(1 + exp(X %*% beta_))

  lpy_c <- sum(log_like(beta * z))
  for (s in seq_len(n_iter)) {
    for (j in sample(2:p)) {
      zp <- z
      zp[j] <- 1 - zp[j]
      lpy_p <- sum(log_like(beta * zp))
      # p(z[j] == 1|y, X, z[-j]) is always numerator
      loj <- (lpy_p - lpy_c) * -1^(zp[j] == 0)
      z[j] <- rbinom(1, 1, 1 / (1 + exp(-loj)))
      # if flipped z is obtained
      if (z[j] == zp[j]) {
        lpy_c <- lpy_p
      }
    }

    beta_p <- t(chol(sigp^2 * diag(p))) %*% rnorm(p) + beta
    # beta_p <- t(mvtnorm::rmvnorm(1, beta, sigp^2 * diag(p)))
    lhr <- sum(log_like(beta_p * z)) - sum(log_like(beta * z)) +
      sum(log_prior(beta_p)) - sum(log_prior(beta))

    if (log(runif(1)) < lhr) {
      beta <- beta_p
      acs <- acs + 1
    }

    BETA[s, ] <- beta
    Z[s, ] <- z
  }
  return(list(acs = acs, BETA = BETA, Gam = Z))
}


n_iter <- 25000
res <- metro_hasting(y, X, 0.07, n_iter)
res$acs / n_iter

thin <- seq(25, n_iter, by = 25)
acf(res$BETA[thin, ])

par(mfrow = c(2, 3))
for (c in seq_len(ncol(X))) {
  plot(res$BETA[thin, c], type = "l")
}
dev.off()

BG <- res$BETA * res$Gam
par(mfrow = c(2, 3))
for (c in seq_len(ncol(X))) {
  plot(BG[thin, c], type = "l")
}
dev.off()


# b --------------------------------------------------

combid <- \(x) sum(c(1, 2, 4, 8, 16, 32) * x)
sort(table(apply(res$Gam, 1, combid)), decreasing = TRUE)
sort(table(apply(res$Gam[thin, ], 1, combid)), decreasing = TRUE)

as.integer(intToBits(63))[1:6]
as.integer(intToBits(59))[1:6]
as.integer(intToBits(57))[1:6]
as.integer(intToBits(61))[1:6]
as.integer(intToBits(27))[1:6]


# c --------------------------------------------------

par(mfrow = c(2, 3))
for (c in seq_len(ncol(X))) {
  plot(density(BG[thin, c]))
}
dev.off()

colnames(X)
apply(BG[thin, ], 2, mean)
apply(res$Gam[thin, ], 2, mean)
