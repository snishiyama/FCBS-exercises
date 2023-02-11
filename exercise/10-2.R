library(coda)
source("sample-code/rmvnorm.R")
dat <- read.table("data/msparrownest.dat", col.names = c("y", "x"))
head(dat)
plot(dat$x, dat$y)

logit <- function(x) log(x / (1 - x))
plot(x <- seq(0.01, 0.99, length.out = 100), logit(x), type = "l")


# b --------------------------------------------------

theta_up <- 0.99
theta_lw <- 0.01
x_up <- 15
x_lw <- 10

beta_up <- (logit(theta_up) - logit(theta_lw)) / (x_up - x_lw)
alp_lw <- logit(theta_lw) - beta_up * x_lw

beta_lw <- (logit(theta_lw) - logit(theta_up)) / (x_up - x_lw)
alp_up <- logit(theta_up) - beta_lw * x_lw

print(c(beta_lw, beta_up))
print(c(alp_lw, alp_up))

m_beta <- (beta_up + beta_lw) / 2
m_alp <- (alp_up + alp_lw) / 2
sig_beta <- (beta_up + m_beta) / 2
sig_alp <- (alp_up + m_alp) / 2

log_prior <- \(beta) dnorm(beta, c(m_alp, m_beta), c(sig_alp, sig_beta), log = TRUE)


# c --------------------------------------------------

metropolis <- \(y, X, n_iter) {
  p <- dim(X)[2]
  var_prop <- 30 * var(y) * solve(t(X) %*% X)
  beta <- rep(0, p)
  acs <- 0
  BETA <- matrix(nrow = n_iter, ncol = p)

  log_like <- \(beta_) y * (X %*% beta_) - log(1 + exp(X %*% beta_))

  for (s in seq_len(n_iter)) {
    beta_p <- t(rmvnorm(1, beta, var_prop))

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

X <- cbind(rep(1, length(dat$x)), dat$x)
n_iter <- 20000
res <- metropolis(dat$y, X, n_iter)

res$acs / n_iter
id_thin10 <- seq(10, n_iter, by = 10)
BETA_thin10 <- res$BETA[id_thin10, ]

png("fig/10-2-c.png", height = 600, width = 1200, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
plot(BETA_thin10[, 2], type = "l", ylab = expression(beta))
acf(BETA_thin10[, 2], main = "")
dev.off()
coda::effectiveSize(BETA_thin10[, 2])
coda::effectiveSize(BETA_thin10[, 1])


# d --------------------------------------------------

png("fig/10-2-d.png", height = 600, width = 1200, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 1))
plot(density(BETA_thin10[, 1], adjust = 2), xlab = "", main = "")
lines(x <- seq(-30, 10, length.out = 1000), dnorm(x, m_alp, sig_alp), lty = "dashed")
title(xlab = expression(alpha))

plot(density(BETA_thin10[, 2], adjust = 2), xlab = "", main = "")
lines(x <- seq(-1, 3, length.out = 1000), dnorm(x, m_beta, sig_beta), lty = "dashed")
title(xlab = expression(beta))
dev.off()


# e --------------------------------------------------

logistic <- \(x) 1 / (1 + exp(-x))
xs <- cbind(rep(1, 100), seq(10, 16, length.out = 100))
theta_pred <- logistic(xs %*% t(BETA_thin10))
theta_q <- apply(theta_pred, 1, quantile, c(0.025, 0.5, 0.975))
head(theta_q)

png("fig/10-2-e.png", height = 600, width = 600, res = 150)
par(mar = c(5, 5, 2, 2))
plot(dat$x, dat$y, xlab = "x", ylab = "y")
lines(xs[, 2], theta_q[1, ], lty = "dashed")
lines(xs[, 2], theta_q[2, ])
lines(xs[, 2], theta_q[3, ], lty = "dashed")
dev.off()

# predictive check
n_sample <- 10000
rb <- BETA_thin10[sample(seq_along(id_thin10), n_sample, replace = TRUE), ]
my_pred <- numeric(n_sample)
for (i in seq_len(n_sample)) {
  rtheta <- c(logistic(X %*% rb[i, ]))
  my_pred[i] <- mean(rbinom(length(rtheta), 1, rtheta))
}
hist(my_pred)
abline(v = mean(dat$y), col = "red", lwd = 3)
box()
