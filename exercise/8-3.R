if (!require(coda)) {
  install.packages("coda")
}
library(coda)

score <- NULL
school_id <- NULL
for (i in 1:8) {
  d <- scan(sprintf("data/school%d.dat", i))
  score <- c(score, d)
  school_id <- c(school_id, rep(i, length(d)))
}

dat <- cbind(school_id, score)


# a --------------------------------------------------

gibbs <- \(data, n_iter = 5000) {
  # prior
  mu0 <- 7
  g20 <- 5
  t20 <- 10
  eta0 <- 2
  s20 <- 15
  nu0 <- 2

  # init
  m <- length(unique(data[, 1])) # n of groups

  ybar <- tapply(data[, 2], data[, 1], mean)
  sv <- tapply(data[, 2], data[, 1], var)
  n <- tapply(data[, 2], data[, 1], length) # n of each group

  theta <- ybar
  sigma2 <- mean(sv)
  mu <- mean(theta)
  tau2 <- var(theta)

  THETA <- matrix(nrow = n_iter, ncol = m)
  SMT <- matrix(nrow = n_iter, ncol = 3)

  for (s in seq_len(n_iter)) {
    # theta
    for (j in seq_len(m)) {
      vtheta <- 1 / (n[j] / sigma2 + 1 / tau2)
      etheta <- vtheta * (n[j] * ybar[j] / sigma2 + mu / tau2)
      theta[j] <- rnorm(1, etheta, sqrt(vtheta))
    }

    # sigma2
    nun <- nu0 + sum(n)
    ss <- nu0 * s20
    for (j in seq_len(m)) {
      ss <- ss + sum((data[data[, 1] == j, 2] - theta[j])^2)
    }
    sigma2 <- 1 / rgamma(1, nun / 2, ss / 2)

    # mu
    vmu <- 1 / (m / tau2 + 1 / g20)
    emu <- vmu * (m * mean(theta) / tau2 + mu0 / g20)
    mu <- rnorm(1, emu, sqrt(vmu))

    # tau2
    etam <- eta0 + m
    ss <- eta0 * t20 + sum((theta - mu)^2)
    tau2 <- 1 / rgamma(1, etam / 2, ss / 2)

    THETA[s, ] <- theta
    SMT[s, ] <- c(sigma2, mu, tau2)
  }

  return(list(THETA = THETA, sigma2 = SMT[, 1], mu = SMT[, 2], tau2 = SMT[, 3]))
}

n_sample <- 5000
post <- gibbs(dat, n_sample)

with(post, lapply(list(sigma2 = sigma2, mu = mu, tau2 = tau2), effectiveSize))


# b --------------------------------------------------

prior <- list(
  sigma2 = 1 / rgamma(n_sample, 1, 15),
  mu = rnorm(n_sample, 7, 5),
  tau2 = 1 / rgamma(n_sample, 1, 10)
)

png("fig/8-3-b.png", width = 800, height = 800, res = 100)
par(mfrow = c(2, 2), mar = c(5, 3, 2, 2))
# sigma2
plot(dens_s2 <- density(post$sigma2, adjust = 2),
  lwd = 3,
  xlab = expression(sigma^2), main = ""
)
lines(density(prior$sigma2, n = max(prior$sigma2) / dens_s2$bw), lwd = 3, col = "gray60")
# mu
plot(density(post$mu, adjust = 2),
  lwd = 3,
  xlab = expression(mu), main = ""
)
lines(density(prior$mu, adjust = 2), lwd = 3, col = "gray60")
# tau2
plot(dens_t2 <- density(post$tau2, adjust = 2),
  lwd = 3,
  xlab = expression(tau^2), main = ""
)
lines(density(prior$tau2, n = max(prior$tau2) / dens_t2$bw), lwd = 3, col = "gray60")
legend("topright",
  legend = c("prior", "posterior"),
  lty = "solid", lwd = c(3, 3), col = c("gray60", "black")
)
dev.off()


# c --------------------------------------------------

png("fig/8-3-c.png", width = 800, height = 800, res = 100)
par(mar = c(5, 3, 2, 3))
plot(dens_r <- with(post, density(tau2 / (tau2 + sigma2), adjust = 2)),
  lwd = 3,
  xlab = expression(R), main = ""
)
lines(with(prior, density(tau2 / (tau2 + sigma2), adjust = 2)),
  lwd = 3,
  col = "gray60"
)
dev.off()


# d --------------------------------------------------

mean(post$THETA[, 6] > post$THETA[, 7])
table(apply(post$THETA, 1, which.min)) / 5000


# e --------------------------------------------------

ybars <- tapply(dat[, 2], dat[, 1], mean)
thetas <- apply(post$THETA, 2, mean)

png("fig/8-3-e.png", width = 600, height = 600, res = 100)
plot(1:8 - 0.1, thetas,
  pch = 19, ylim = c(6, 11),
  xlab = "School", ylab = "Mean homework time per week"
)
points(1:8 + 0.1, ybars, pch = 24)
abline(h = mean(dat[, 2]), lty = "dashed")
abline(h = mean(post$mu))
legend("topleft",
  legend = c("posterior", "data"),
  pch = c(19, 24), lty = c("solid", "dashed")
)
dev.off()
