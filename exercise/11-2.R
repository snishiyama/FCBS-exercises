library(mvtnorm)
library(monomvn)
dat <- read.table("data/pdensity.dat", header = TRUE)
head(dat)

dat$density2 <- dat$density^2
head(dat)


# a --------------------------------------------------

n_groups <- length(unique(dat$plot))
Beta_ols <- matrix(ncol = 3, nrow = n_groups)
resid_ols <- numeric(n_groups)
for (p in unique(dat$plot)) {
  dat_p <- dat[dat$plot == p, ]
  lm_fit <- lm(yield ~ density + density2, data = dat_p)
  Beta_ols[p, ] <- lm_fit$coefficients
  resid_ols[p] <- (summary(lm_fit)$sigma)^2
}

png("fig/11-2-a.png", res = 100)
# par(mar = c(4.5, 4.5, 1, 1))
plot(dat$density, dat$yield,
  xlab = "Density", ylab = "Yield",
  main = "OLS regression"
)
for (i in seq_len(nrow(Beta_ols))) {
  regline <- \(x) Beta_ols[i, 1] + Beta_ols[i, 2] * x + Beta_ols[i, 3] * x^2
  lines(x <- seq(1, 9, length.out = 50), regline(x))
}
dev.off()

Beta_ols
(theta_hat <- apply(Beta_ols, 2, mean))
# following two lines are equivalent to cov(Beta_ols)
Beta_olsc <- t(t(Beta_ols) - theta_hat)
(Sigma_hat <- t(Beta_olsc) %*% Beta_olsc / 9)
(sig2_hat <- mean(resid_ols))


# b, c --------------------------------------------------

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

# y is n (dim) vector
# X is n x p matrix (p is # of iv, n is # of observation)
# Beta is p x g matrix (g is # of groups)
# `groups` is n vector of group ids
sample_sig2 <- \(y, X, groups, Beta, nu0, sig20) {
  n <- length(y)
  SSR <- 0
  for (g in unique(groups)) {
    yg <- y[groups == g]
    Xg <- X[groups == g, ]
    betag <- Beta[, unique(groups) == g]
    SSR <- SSR + sum((yg - Xg %*% betag)^2)
  }
  return(1 / rgamma(1, (nu0 + n) / 2, (nu0 * sig20 + SSR) / 2))
}

sample_Beta <- \(y, X, groups, theta, Sigma, sig2) {
  Beta <- NULL
  iSigma <- solve(Sigma)
  for (g in unique(groups)) {
    Xg <- X[groups == g, ]
    yg <- y[groups == g]
    vBeta <- solve(iSigma + t(Xg) %*% Xg / sig2)
    eBeta <- vBeta %*% (iSigma %*% theta + t(Xg) %*% yg / sig2)
    beta <- mvtnorm::rmvnorm(1, eBeta, vBeta)
    Beta <- rbind(Beta, beta)
  }
  return(t(Beta))
}

gibbs <- \(y, X, groups, mu0, Lambda0, eta0, S0, nu0, sig20, n_iter) {
  p <- ncol(X)
  n_group <- length(unique(groups))

  post_samples <- list(
    THETA = matrix(nrow = n_iter, ncol = p),
    SIGMA = matrix(nrow = n_iter, ncol = p * p),
    BETA = array(dim = c(n_iter, p, n_group)), # 3d array
    SIG2 = numeric(n_iter)
  )

  # init
  theta <- mu0
  Sigma <- S0
  sig2 <- sig20
  Beta <- matrix(nrow = p, ncol = n_group)
  for (i in seq_len(n_group)) {
    g <- unique(groups)[i]
    Xg <- X[groups == g, ]
    yg <- y[groups == g]
    Beta[, i] <- solve(t(Xg) %*% Xg) %*% t(Xg) %*% yg
  }

  # main gibbs sampling
  stime <- Sys.time()
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) {
      dt <- round(difftime(Sys.time(), stime, units = "secs"), 2)
      cat("Running", i, "th loop. Elapsed", dt, "sec\n")
    }
    Beta <- sample_Beta(y, X, groups, theta, Sigma, sig2)
    theta <- sample_theta(Beta, Sigma, mu0, Lambda0)
    Sigma <- sample_Sigma(Beta, theta, eta0, S0)
    sig2 <- sample_sig2(y, X, groups, Beta, nu0, sig20)

    post_samples$BETA[i, , ] <- Beta
    post_samples$THETA[i, ] <- t(theta)
    post_samples$SIGMA[i, ] <- t(c(Sigma))
    post_samples$SIG2[i] <- sig2
  }
  return(post_samples)
}

X <- as.matrix(dat[, c("density", "density2")])
X <- cbind(rep(1, nrow(X)), X)
colnames(X)[1] <- "intercept"
head(X)

# about 30 s take
n_iter <- 20000
res <- gibbs(
  dat$yield, X, dat$plot,
  theta_hat, Sigma_hat, # prior for theta
  4, Sigma_hat, # for Sigma
  2, sig2_hat, # for sig2
  n_iter
)

# par(mfrow = c(2, 2))
# for (i in 1:3) {
#   plot(res$THETA[, i], type = "l")
# }
# plot(res$SIG2, type = "l")

lag <- 10
thinning <- seq(lag, n_iter, by = lag)
post_thin <- list(
  BETA = res$BETA[thinning, , ],
  THETA = res$THETA[thinning, ],
  SIGMA = res$SIGMA[thinning, ],
  SIG2 = res$SIG2[thinning]
)

# effective size
coda::effectiveSize(post_thin$THETA)
coda::effectiveSize(post_thin$SIGMA)
coda::effectiveSize(post_thin$SIG2)
apply(post_thin$BETA, 3, coda::effectiveSize)

# autocorrelation
acf(post_thin$THETA)
acf(post_thin$SIGMA[, 1:5])
acf(post_thin$SIGMA[, 6:9])
acf(post_thin$SIG2)
acf(post_thin$BETA[, , 1])
dim(post_thin$BETA)


(Beta_bayes <- t(apply(post_thin$BETA, c(2, 3), mean)))

png("fig/11-2-c.png", res = 100)
plot(dat$density, dat$yield,
  xlab = "Density", ylab = "Yield",
  main = "Bayesian regression"
)
for (i in seq_len(nrow(Beta_bayes))) {
  reg_bayes <- \(x) Beta_bayes[i, 1] + Beta_bayes[i, 2] * x + Beta_bayes[i, 3] * x^2
  lines(x <- seq(1, 9, length.out = 50), reg_bayes(x))
}
dev.off()


# d --------------------------------------------------

THETA_prior <- mvtnorm::rmvnorm(n_iter, theta_hat, Sigma_hat)
png("fig/11-2-d_theta.png", height = 800, width = 800, res = 120)
par(mfrow = c(2, 2))
for (i in 1:3) {
  plot(density(post_thin$THETA[, i], adjust = 2), main = bquote(theta[.(i)]))
  lines(density(THETA_prior[, i], adjust = 2), lty = "dashed")
}
dev.off()

SIGMA_prior <- matrix(nrow = n_iter, ncol = 9)
for (i in seq_len(n_iter)) {
  SIGMA_prior[i, ] <- t(c(solve(monomvn::rwish(4, solve(Sigma_hat)))))
}
png("fig/11-2-d_sigma.png", height = 1200, width = 1200, res = 120)
par(mfrow = c(3, 3))
for (i in seq_len(ncol(post_thin$SIGMA))) {
  i_ <- floor((i - 1) / 3) + 1
  j <- floor((i - 1) %% 3) + 1
  sigma <- post_thin$SIGMA[, i]
  plot(dens_post <- density(sigma, adjust = 2),
    xlim = quantile(sigma, c(0.001, 0.999)),
    main = bquote(Sigma[.(i_) * "," * .(j)])
  )
  sig_prior <- SIGMA_prior[, i]
  # for compute # of bins. abs() is required for negative values
  edge <- max(abs(range(sig_prior)))
  n_bin <- as.integer(edge / dens_post$bw)
  lines(density(sig_prior, bw = dens_post$bw, n = n_bin), lty = "dashed")
}
dev.off()


# e --------------------------------------------------

xmax <- \(theta) round(theta[2] / (-2 * theta[3]))

table(apply(post_thin$THETA, 1, xmax)) / length(thinning)

pred <- numeric(nrow(post_thin$THETA))
for (i in seq_len(nrow(post_thin$THETA))) {
  theta <- post_thin$THETA[i, ]
  p <- length(theta)
  sigma <- matrix(post_thin$SIGMA[i, ], p, p)
  beta <- t(mvtnorm::rmvnorm(1, theta, sigma))
  x <- c(1, 6, 6^2)
  pred[i] <- rnorm(1, t(x) %*% beta, sqrt(post_thin$SIG2[i]))
}

quantile(pred, c(0.025, 0.5, 0.975))
plot(density(pred, adjust = 2))
