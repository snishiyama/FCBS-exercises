library(mvtnorm)
library(monomvn)
library(coda)
dat <- read.table("data/mathstandard2.dat", header = TRUE)
head(dat)

logistic <- \(x) 1 / (1 + exp(-x))


# b --------------------------------------------------

counties <- unique(dat$county)

Beta_ml <- matrix(nrow = length(counties), ncol = 2)
for (i in seq_len(length(counties))) {
  c <- counties[i]
  dat_c <- dat[dat$county == c, ]
  yj <- dat_c$metstandard
  xj <- dat_c$percentms
  Beta_ml[i, ] <- glm(yj ~ xj, family = binomial)$coef
}
Beta_ml

glmline <- \(beta, x) logistic(beta[1] + beta[2] * x)

png("fig/11-4-b.png", width = 1400, height = 1200, res = 100)
par(mfrow = c(6, 7), mar = c(2.5, 2.5, 3, 1))
for (i in seq_len(length(counties))) {
  c <- counties[i]
  dat_c <- dat[dat$county == c, ]
  yj <- dat_c$metstandard
  xj <- dat_c$percentms
  plot(xj, yj, main = c, ylim = c(0, 1))

  beta <- Beta_ml[i, ]
  if (any(is.na(beta))) {
    next
  }
  lines(x <- seq(min(xj), max(xj), length.out = 50), glmline(beta, x))
}
dev.off()

(n_school <- sapply(counties, \(c) sum(dat$county == c)))
(theta_hat <- apply(Beta_ml[n_school > 9, ], 2, mean))
(Sigma_hat <- var(Beta_ml[n_school > 9, ]))


# c --------------------------------------------------

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
sample_Beta <- \(y, X, groups, theta, Sigma, Beta, adj_prop) {
  log_prior <- \(beta) mvtnorm::dmvnorm(t(beta), theta, Sigma, log = TRUE)
  log_like <- \(beta, y_, X_) y_ * (X_ %*% beta) - log(1 + exp(X_ %*% beta))

  Beta_new <- NULL
  Acs <- NULL
  ugroups <- unique(groups)
  # metropolis for each group
  for (g in ugroups) {
    Xg <- X[groups == g, ]
    yg <- y[groups == g]
    betag <- Beta[, ugroups == g]

    beta_p <- t(rmvnorm(1, betag, adj_prop * Sigma))

    lhr <- sum(log_like(beta_p, yg, Xg)) - sum(log_like(betag, yg, Xg)) +
      log_prior(beta_p) - log_prior(betag)

    acs <- 0
    beta <- betag
    if (log(runif(1)) < lhr) {
      beta <- beta_p
      acs <- 1
    }

    Beta_new <- rbind(Beta_new, c(beta))
    Acs <- c(Acs, acs)
  }
  return(list(Beta = t(Beta_new), Acs = Acs))
}

mhmcmc <- \(y, X, groups, mu0, Lambda0, eta0, S0, adj_prop, n_iter) {
  p <- ncol(X)
  n_group <- length(unique(groups))

  post_samples <- list(
    THETA = matrix(nrow = n_iter, ncol = p),
    SIGMA = matrix(nrow = n_iter, ncol = p * p),
    BETA = array(dim = c(n_iter, p, n_group)), # 3d array
    ACS = numeric(n_group)
  )

  # init
  theta <- mu0
  Sigma <- S0
  Beta <- matrix(theta, nrow = p, ncol = n_group)

  # main gibbs sampling
  start <- Sys.time()
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) {
      diff <- round(difftime(Sys.time(), start, unit = "secs"), 2)
      cat("Running", i, "th loop. Elapsed", diff, "sec\n")
    }
    theta <- sample_theta(Beta, Sigma, mu0, Lambda0)
    Sigma <- sample_Sigma(Beta, theta, eta0, S0)
    beta_acs <- sample_Beta(y, X, groups, theta, Sigma, Beta, adj_prop)
    Beta <- beta_acs$Beta

    post_samples$BETA[i, , ] <- Beta
    post_samples$THETA[i, ] <- t(theta)
    post_samples$SIGMA[i, ] <- t(c(Sigma))
    post_samples$ACS <- post_samples$ACS + beta_acs$Acs
  }
  post_samples$ACS <- post_samples$ACS / n_iter
  return(post_samples)
}

X <- as.matrix(dat[, "percentms", drop = FALSE])
X <- cbind(rep(1, length(X)), X)
colnames(X)[1] <- "intercept"
head(X)

fnm_mcmc <- "data/11-4_mcmcsample.rds"
if (file.exists(fnm_mcmc)) {
  res_thin <- readRDS(fnm_mcmc)
} else {
  # About 12 mins in my environment
  n_iter <- 100000
  res <- mhmcmc(
    dat$metstandard, X, dat$county,
    theta_hat, Sigma_hat, # prior for theta
    4, Sigma_hat, # for Sigma
    1.2,
    n_iter
  )

  lag <- 40
  thinning <- seq(lag, n_iter, by = lag)
  res_thin <- list(
    THETA = res$THETA[thinning, ],
    SIGMA = res$SIGMA[thinning, ],
    BETA = res$BETA[thinning, , ],
    ACS = res$ACS
  )
  saveRDS(res_thin, fnm_mcmc)
}

summary(res_thin$ACS)

# d --------------------------------------------------

# diagnosis
acf(res_thin$THETA)
acf(res_thin$SIGMA[, c(1, 2, 4)])
acf(res_thin$BETA[, , 1])

coda::effectiveSize(res_thin$THETA)
coda::effectiveSize(res_thin$SIGMA)
apply(res_thin$BETA, 3, coda::effectiveSize)

png("fig/11-4-d.png", height = 800, width = 1200, res = 100)
par(mfrow = c(2, 3))
plot(density(res_thin$THETA[, 1], adjust = 2), main = expression(theta[1]))
plot(density(res_thin$THETA[, 2], adjust = 2), main = expression(theta[2]))
plot.new()
plot(density(res_thin$SIGMA[, 1], adjust = 2), main = expression(Sigma[1 * "," * 1]))
plot(density(res_thin$SIGMA[, 2], adjust = 2), main = expression(Sigma[1 * "," * 2]))
plot(density(res_thin$SIGMA[, 4], adjust = 2), main = expression(Sigma[2 * "," * 2]))
dev.off()

# e --------------------------------------------------

(Beta_bayes <- t(apply(res_thin$BETA, c(2, 3), mean)))

png("fig/11-4-e.png", width = 1400, height = 1200, res = 100)
par(mfrow = c(6, 7), mar = c(2.5, 2.5, 3, 1))
for (i in seq_len(length(counties))) {
  c <- counties[i]
  dat_c <- dat[dat$county == c, ]
  yj <- dat_c$metstandard
  xj <- dat_c$percentms
  plot(xj, yj, main = c, ylim = c(0, 1))

  lines(x <- seq(10, 110, length.out = 50), glmline(Beta_bayes[i, ], x))
  beta <- Beta_ml[i, ]
  if (any(is.na(beta))) {
    next
  }
  lines(x, glmline(beta, x), lty = "dashed")
}
plot.new()
legend("center",
  legend = c("Bayes", "ML"), lty = c("solid", "dashed"),
  cex = 2.3, box.lty = "blank"
)
dev.off()

# f --------------------------------------------------

theta_prior <- mvtnorm::rmvnorm(10000, theta_hat, Sigma_hat)
Sigma_prior <- matrix(nrow = 10000, ncol = 4)
for (i in seq_len(10000)) {
  Sigma_prior[i, ] <- t(c(solve(monomvn::rwish(4, solve(Sigma_hat)))))
}

png("fig/11-4-f.png", height = 800, width = 1200, res = 120)
par(mfrow = c(2, 3))
for (i in 1:2) {
  plot(density(res_thin$THETA[, i], adjust = 2), main = bquote(theta[.(i)]))
  lines(density(theta_prior[, i], adjust = 2), lty = "dashed")
  abline(v = theta_hat[i], col = "blue")
}
plot.new()
for (i in 1:4) {
  if (i == 3) next
  i_ <- trunc((i - 1) / 2) + 1
  j <- (i - 1) %% 2 + 1
  plot(dens <- density(res_thin$SIGMA[, i], adjust = 2),
    main = bquote(Sigma[.(i_) * "," * .(j)])
  )
  lines(density(sp <- Sigma_prior[, i], n = max(abs(sp)) / dens$bw), lty = "dashed")
  abline(v = Sigma_hat[i_, j], col = "blue")
}
dev.off()
