library(mvtnorm)
library(monomvn)
dat <- read.table("data/pdensity.dat", header = TRUE)
head(dat)

dat$density2 <- dat$density^2
head(dat)

n_groups <- length(unique(dat$plot))
Beta_ols <- matrix(ncol = 3, nrow = n_groups)
resid_ols <- numeric(n_groups)
for (p in unique(dat$plot)) {
  dat_p <- dat[dat$plot == p, ]
  lm_fit <- lm(yield ~ density + density2, data = dat_p)
  Beta_ols[p, ] <- lm_fit$coefficients
  resid_ols[p] <- (summary(lm_fit)$sigma)^2
}
theta_hat <- apply(Beta_ols, 2, mean)
Sigma_hat <- cov(Beta_ols)


# a --------------------------------------------------

# Same as the full conditional dist. in Section 8.5.
# sig2_ is a vector of groups' sig2, m (dim) vector (m is # of groups)
sample_sig20 <- \(sig2_, nu0, a, b) {
  m <- length(sig2_)
  return(rgamma(1, a + m * nu0 / 2, b + nu0 * sum(1 / sig2_) / 2))
}


# b --------------------------------------------------

# Same as the full conditional dist. in Section 8.5,
# except that theta_j is replaced with t(beta_j) * x_{i, j}
# y is n (dim) vector
# X is n x p matrix (p is # of iv, n is # of observation)
# Beta is p x g matrix (g is # of groups)
# `groups` is n vector of group ids
sample_sig2_ <- \(y, X, groups, Beta, nu0, sig20) {
  sig2_ <- NULL
  for (g in unique(groups)) {
    yg <- y[groups == g]
    Xg <- X[groups == g, ]
    betag <- Beta[, unique(groups) == g]

    SSRg <- sum((yg - Xg %*% betag)^2)
    sig2g <- 1 / rgamma(1, (nu0 + length(yg)) / 2, (nu0 * sig20 + SSRg) / 2)
    sig2_ <- c(sig2_, sig2g)
  }
  return(sig2_)
}


# c --------------------------------------------------

# Same as the full conditional dist. in Section 11.2,
# except that sigma^2 is replaced with sigma^2_j
sample_Beta <- \(y, X, groups, theta, Sigma, sig2_) {
  Beta <- NULL
  iSigma <- solve(Sigma)
  ugroups <- unique(groups)
  for (g in ugroups) {
    Xg <- X[groups == g, ]
    yg <- y[groups == g]
    sig2g <- sig2_[ugroups == g]
    vBeta <- solve(iSigma + t(Xg) %*% Xg / sig2g)
    eBeta <- vBeta %*% (iSigma %*% theta + t(Xg) %*% yg / sig2g)
    beta <- mvtnorm::rmvnorm(1, eBeta, vBeta)
    Beta <- rbind(Beta, beta)
  }
  return(t(Beta))
}


# d --------------------------------------------------

# sig2_ is a vector of groups' sig2, m (dim) vector (m is # of groups)
# lr_nu0 <- \(nu0p, nu0s, sig20, sig2_) {
#   m <- length(sig2_)
#   sl_isig2j <- sum(log(1 / sig2_))
#   sum_isig2j <- sum(1 / sig2_)

#   return(
#     m * (0.5 * (nu0p * log(nu0p * sig20 / 2) - nu0s * log(nu0s * sig20 / 2)) -
#       (lgamma(0.5 * nu0p) - lgamma(0.5 * nu0s))) +
#       0.5 * (nu0p - nu0s) * sl_isig2j -
#       0.5 * (nu0p - nu0s) * sig20 * sum_isig2j
#   )
# }

lr_nu0 <- \(nu0p, nu0s, sig20, sig2_) {
  return(
    sum(dgamma(1 / sig2_, nu0p / 2, nu0p * sig20 / 2, log = TRUE)) -
      sum(dgamma(1 / sig2_, nu0s / 2, nu0s * sig20 / 2, log = TRUE))
  )
}

# metropolis-hastings for nu0
sample_nu0 <- \(nu0s, sig20, sig2_, delta, nu0max = 100) {
  nu0p <- sample(seq(nu0s - delta, nu0s + delta, by = 1), 1)
  if (nu0p > nu0max) {
    nu0p <- (2 * nu0max + 1) - nu0p
  } else if (nu0p < 1) {
    nu0p <- 1 - nu0p
  }
  lr <- lr_nu0(nu0p, nu0s, sig20, sig2_)
  if (log(runif(1)) < lr) {
    return(list(nu0 = nu0p, acs = 1))
  }
  return(list(nu0 = nu0s, acs = 0))
}


# e --------------------------------------------------

# Beta is p x g matrix (p is # of iv, g is # of groups)
# Sigma is p x p matrix
sample_theta <- \(Beta, Sigma, mu0, Lambda0) {
  n_group <- ncol(Beta)
  iSigma <- solve(Sigma)
  iLambda0 <- solve(Lambda0)

  mBeta <- apply(Beta, 1, mean)
  Lambdam <- solve(iLambda0 + n_group * iSigma)
  mum <- Lambdam %*% (iLambda0 %*% mu0 + n_group * iSigma %*% mBeta)
  return(t(mvtnorm::rmvnorm(1, mum, Lambdam)))
}

# Beta is p x g matrix (p is # of iv, g is # of groups)
# theta is p (dim) vector
sample_Sigma <- \(Beta, theta, eta0, S0) {
  g <- ncol(Beta)
  Sn <- S0 + (Beta - c(theta)) %*% t(Beta - c(theta))
  return(solve(monomvn::rwish(eta0 + g, solve(Sn))))
}

mhmcmc <- \(y, X, groups, mu0, Lambda0, eta0, S0, a, b, delta, nu0max, n_iter) {
  p <- ncol(X)
  n_group <- length(unique(groups))

  post_samples <- list(
    THETA = matrix(nrow = n_iter, ncol = p),
    SIGMA = matrix(nrow = n_iter, ncol = p * p),
    SIG20 = numeric(n_iter),
    NU0 = numeric(n_iter),
    BETA = array(dim = c(n_iter, p, n_group)), # 3d array
    SIG2_ = matrix(nrow = n_iter, ncol = n_group),
    ACS = 0
  )


  # init
  theta <- mu0
  Sigma <- S0
  sig20 <- (summary(lm(y ~ -1 + X))$sigma)^2
  nu0 <- sample(1:nu0max, 1)
  Beta <- matrix(nrow = p, ncol = n_group)
  sig2_ <- numeric(n_group)
  for (i in seq_len(n_group)) {
    g <- unique(groups)[i]
    Xg <- X[groups == g, ]
    yg <- y[groups == g]
    Beta[, i] <- solve(t(Xg) %*% Xg) %*% t(Xg) %*% yg
    sig2_[i] <- (summary(lm(yg ~ Xg))$sigma)^2
  }

  # main sampling
  start <- Sys.time()
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) {
      diff <- round(difftime(Sys.time(), start, unit = "secs"), 2)
      cat("Running", i, "th loop. Elapsed", diff, "sec\n")
    }
    # hyper
    theta <- sample_theta(Beta, Sigma, mu0, Lambda0)
    Sigma <- sample_Sigma(Beta, theta, eta0, S0)
    sig20 <- sample_sig20(sig2_, nu0, a, b)
    nu0_acs <- sample_nu0(nu0, sig20, sig2_, delta, nu0max)
    nu0 <- nu0_acs$nu0

    # groups
    Beta <- sample_Beta(y, X, groups, theta, Sigma, sig2_)
    sig2_ <- sample_sig2_(y, X, groups, Beta, nu0, sig20)

    post_samples$THETA[i, ] <- t(theta)
    post_samples$SIGMA[i, ] <- t(c(Sigma))
    post_samples$SIG20[i] <- sig20
    post_samples$NU0[i] <- nu0
    post_samples$BETA[i, , ] <- Beta
    post_samples$SIG2_[i, ] <- sig2_
    post_samples$ACS <- post_samples$ACS + nu0_acs$acs
  }
  post_samples$ACS <- post_samples$ACS / n_iter
  return(post_samples)
}

X <- as.matrix(dat[, c("density", "density2")])
X <- cbind(rep(1, nrow(X)), X)
colnames(X)[1] <- "intercept"
head(X)

n_iter <- 20000
res <- mhmcmc(
  dat$yield, X, dat$plot,
  theta_hat, Sigma_hat, # prior for theta
  4, Sigma_hat, # for Sigma
  2, 2, # for sig20
  80, 100, # for nu0
  n_iter
)
res$ACS

lag <- 10
thinning <- seq(lag, n_iter, by = lag)
res_thin <- list(
  THETA = res$THETA[thinning, ],
  SIGMA = res$SIGMA[thinning, ],
  SIG20 = res$SIG20[thinning],
  NU0 = res$NU0[thinning],
  BETA = res$BETA[thinning, , ],
  SIG2_ = res$SIG2_[thinning, ]
)

acf(res_thin$NU0)
coda::effectiveSize(res_thin$NU0)

plot(res_thin$SIG20, type = "l")
plot(res_thin$NU0, type = "l")


acf(res_thin$THETA)
acf(res_thin$SIGMA[, 1:5])
acf(res_thin$SIGMA[, 6:9])
acf(res_thin$SIG20)
acf(res_thin$NU0)

plot(res_thin$BETA[, 1, 1], type = "l")
plot(res_thin$SIG2_[, 1], type = "l")


# f --------------------------------------------------

cnt_NU0 <- (table(c(res_thin$NU0, 1:100)) - rep(1, 100)) / length(thinning)
sig2j <- with(res_thin, 1 / rgamma(n_iter / lag, NU0 / 2, NU0 * SIG20 / 2))
summary(sig2j)
summary(res_thin$SIG20)
quantile(sig2j, c(0.025, 0.5, 0.975))

png("fig/11-3-f_nu0.png", width = 1200, height = 600, res = 100)
par(mfrow = c(1, 2))
# prior and posterior nu0
plot(1:100 + 0.1, cnt_NU0,
  type = "h", ylim = range(cnt_NU0),
  main = expression(paste("Prior and posterior dist. of ", nu[0])),
  xlab = expression(nu[0]), ylab = "Density"
)
lines(1:100 - 0.1, rep(1 / 100, 100), type = "h", col = "blue")
axis(2, at = seq(0, max(cnt_NU0), by = 0.005))
legend("topright",
  legend = c("Prior", "Posterior"),
  lty = "solid", col = c("blue", "black")
)
# sig20 and predictive sig2j
plot(density(sig2j, adjust = 2),
  main = expression(paste("Posterior dist. of ", sigma[0]^2, " and predictive dist. of ", sigma[j]^2)),
  xlab = "",
  xlim = c(0, quantile(sig2j, 0.99)), ylim = c(0, 3)
)
lines(density(res_thin$SIG20, adjust = 2), lty = "dashed")
legend("topright",
  legend = c(
    expression(p(tilde(sigma)[j]^2 ~ "|" ~ Y)),
    expression(p(sigma[0]^2 ~ "|" ~ Y))
  ),
  lty = c("solid", "dashed")
)
dev.off()
