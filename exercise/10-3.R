dat <- read.table("data/tplant.dat", col.names = c("height", "time", "pH"))
str(dat)

# a --------------------------------------------------
fit_ols <- lm(height ~ time + pH, data = dat)
summary(fit_ols)


# b --------------------------------------------------

par(mfrow = c(2, 2))
plot(fit_ols)
dev.off()

resi <- summary(fit_ols)$residual
time1 <- seq(1, 19, by = 2)
time2 <- seq(2, 20, by = 2)
plot(resi[time1], resi[time2])
plot(density(resi))

# c
paired_vcov <- \(rho, n) {
  a <- matrix(seq_len(n), ncol = 2, byrow = TRUE)
  v <- matrix(0, n, n)
  for (i in seq_len(n / 2)) {
    ids <- a[i, ]
    v[ids[1], ids[2]] <- rho
  }
  v <- v + t(v)
  diag(v) <- 1
  v
}

mcmc <- \(y, X, delta, n_iter = 10000) {
  n <- length(y)
  p <- ncol(X)
  # prior
  beta_0 <- rep(0, p)
  Sig_0 <- diag(1000, p, p)
  nu0 <- 1
  sig2_0 <- 1

  iSig_0 <- solve(Sig_0)
  XX <- t(X) %*% X

  # init
  Beta <- matrix(nrow = n_iter, ncol = ncol(X))
  colnames(Beta) <- colnames(X)
  sig2 <- numeric(n_iter)
  Rho <- numeric(n_iter)
  acs <- 0

  B <- beta_0
  s2 <- var(y)
  rho <- 0.5
  # main loop
  for (i in seq_len(n_iter)) {
    Crho <- paired_vcov(rho, n)
    iCrho <- solve(Crho)

    # Beta
    V <- solve(iSig_0 + t(X) %*% iCrho %*% X / s2)
    E <- V %*% (iSig_0 %*% beta_0 + t(X) %*% iCrho %*% y / s2)
    B <- t(mvtnorm::rmvnorm(1, E, V))

    # sig2
    dd <- y - X %*% B
    SSRrho <- t(dd) %*% iCrho %*% (dd)
    s2 <- 1 / rgamma(1, (nu0 + n) / 2, (nu0 * sig2_0 + SSRrho) / 2)

    # rho
    rho_p <- abs(runif(1, rho - delta, rho + delta))
    if (rho_p > 1) {
      rho_p <- 2 - rho_p
    }
    Crho_p <- paired_vcov(rho_p, n)
    lr <- mvtnorm::dmvnorm(y, X %*% B, s2 * Crho_p, log = TRUE) -
      mvtnorm::dmvnorm(y, X %*% B, s2 * Crho, log = TRUE)
    if (log(runif(1)) < lr) {
      rho <- rho_p
      acs <- acs + 1
    }

    Beta[i, ] <- B
    sig2[i] <- s2
    Rho[i] <- rho
  }

  return(list(beta = Beta, sig2 = sig2, rho = Rho, acs = acs))
}

t(1:5)
diag(1, 2, 2)

X <- as.matrix(dat[, c("time", "pH")])
X <- cbind(rep(1, nrow(X)), X)
head(X)

res <- mcmc(dat$height, X, 0.3, 10000)

res$acs
thinning <- seq(10, 10000, by = 10)
acf(res$rho[thinning])
plot(res$rho[thinning], type = "l")
plot(density(res$rho[thinning]))

plot(density(res$beta[, 1]))
plot(density(res$beta[, 2]))
plot(density(res$beta[, 3]))

plot(density(res$sig2))

Ebeta <- apply(res$beta[thinning, ], 2, mean)
Esig2 <- mean(res$sig2[thinning])
Erho <- mean(res$rho[thinning])

t(Ebeta)

resi_bayes <- dat$height - X %*% Ebeta
plot(resi_bayes[time1], resi_bayes[time2])

apply(res$beta, 2, quantile, c(0.025, 0.5, 0.975))
quantile(res$sig2, c(0.025, 0.5, 0.975))
