dat <- read.table("data/azdiabetes.dat", header = TRUE)
head(dat)

X <- (\() {
  X <- as.matrix(dat[, c("npreg", "bp", "skin", "bmi", "ped", "age")])
  X <- cbind(rep(1, nrow(X)), X)
  colnames(X)[1] <- "intercept"
  return(X)
})()
y <- dat$glu

regmc <- \(y, X, g, nu0, s20, S) {
  n <- nrow(X)
  p <- ncol(X)
  Hg <- (g / (g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
  SSRg <- t(y) %*% (diag(1, nrow = n) - Hg) %*% y
  # sigma
  s2 <- 1 / rgamma(S, (nu0 + n) / 2, (nu0 * s20 + SSRg) / 2)
  # beta
  Vb <- g / (g + 1) * solve(t(X) %*% X)
  Eb <- Vb %*% t(X) %*% y
  E <- matrix(rnorm(S * p, 0, sqrt(s2)), S, p)
  beta <- t(t(E %*% chol(Vb)) + c(Eb))

  return(list(beta = beta, sigma2 = s2))
}


# a --------------------------------------------------

fit <- regmc(y, X, length(y), 2, 1, 1000)
head(fit$beta)

apply(fit$beta, 2, quantile, probs = c(0.025, 0.5, 0.975))
quantile(fit$sigma2, c(0.025, 0.5, 0.975))

par(mfrow = c(2, 4))
for (nm in colnames(fit$beta)) {
  plot(density(fit$beta[, nm], adjust = 2), xlab = nm, main = "")
  abline(v = 0, lty = "dashed")
}
dev.off()


# b --------------------------------------------------

lpy_X <- \(y, X, g, nu0, s20) {
  n <- nrow(X)
  p <- ncol(X)
  if (p == 0) {
    Hg <- 0
  } else if (p > 0) {
    Hg <- (g / (g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
  }
  SSRg <- t(y) %*% (diag(1, nrow = n) - Hg) %*% y
  lp1 <- -0.5 * n * log(pi)
  lp2 <- lgamma((nu0 + n) / 2) - lgamma(nu0 / 2)
  lp3 <- -0.5 * log(1 + g) * p
  lp4 <- 0.5 * (log(nu0 * s20) * nu0 - log(nu0 * s20 + SSRg) * (nu0 + n))
  return(lp1 + lp2 + lp3 + lp4)
}

# model averaging in 9.3.2
gibbs <- \(y, X, g, nu0, n_iter) {
  z <- rep(1, ncol(X))
  s20 <- (summary(lm(y ~ -1 + X))$sigma)^2
  lpy_c <- lpy_X(y, X, g, nu0, s20)
  Beta <- Z <- matrix(nrow = n_iter, ncol = ncol(X))
  for (i in seq_len(n_iter)) {
    if (i %% 100 == 1) {
      cat("running", i, "st iter\n")
    }
    # z
    for (j in sample(seq_len(ncol(X)))) {
      zp <- z
      zp[j] <- 1 - zp[j]
      Xz <- X[, zp == 1, drop = FALSE]
      if (ncol(Xz) > 0) {
        s20 <- (summary(lm(y ~ -1 + Xz))$sigma)^2
      } else {
        s20 <- mean(y^2)
      }
      lpy_p <- lpy_X(y, Xz, g, nu0, s20)
      # p(z[j] == 1|y, X, z[-j]) is always numerator
      loj <- (lpy_p - lpy_c) * -1^(zp[j] == 0)
      z[j] <- rbinom(1, 1, 1 / (1 + exp(-loj)))
      # if flipped z is obtained
      if (z[j] == zp[j]) {
        lpy_c <- lpy_p
      }
    }
    # beta
    beta <- z
    if (sum(beta) > 0) {
      Xz <- X[, z == 1, drop = FALSE]
      if (ncol(Xz) > 0) {
        s20 <- (summary(lm(y ~ -1 + Xz))$sigma)^2
      } else {
        s20 <- mean(y^2)
      }
      beta[z == 1] <- regmc(y, Xz, g, nu0, s20, 1)$beta
    }
    Z[i, ] <- z
    Beta[i, ] <- beta
  }
  return(list(Z = Z, Beta = Beta))
}

res <- gibbs(y, X, length(y), 1, 10000)

colnames(res$Beta) <- colnames(fit$beta)
colnames(res$Z) <- colnames(fit$beta)
apply(res$Z, 2, mean)
apply(res$Beta, 2, quantile, prob = c(0.025, 0.5, 0.975))

par(mfrow = c(2, 4))
for (nm in colnames(res$Beta)) {
  plot(density(res$Beta[, nm], adjust = 2), xlab = nm, main = "")
  abline(v = 0, lty = "dashed")
}
dev.off()

# model selection in 9.3.1
lp_all_models <- \(y, X, g, nu0) {
  p <- ncol(X) - 1
  zs <- matrix(nrow = 2^p, ncol = ncol(X))
  lps <- numeric(2^p)
  for (i in 0:2^p - 1) {
    zs[i + 1, ] <- c(1, as.integer(intToBits(i))[1:p])
  }
  for (i in seq_len(nrow(zs))) {
    z <- zs[i, ]
    Xz <- X[, z == 1, drop = FALSE]
    s20 <- (summary(lm(y ~ -1 + Xz))$sigma)^2
    lps[i] <- lpy_X(y, Xz, g, nu0, s20)
  }
  return(list(Z = zs, lp = lps))
}

model_comp <- lp_all_models(y, X, length(y), 1)
lp_max <- max(model_comp$lp)
p_model <- exp(model_comp$lp - lp_max) / sum(exp(model_comp$lp - lp_max))

apply(model_comp$Z, 2, \(x) sum(x * p_model))

model_comp$Z[which.max(p_model), ]
