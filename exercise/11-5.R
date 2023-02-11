library(coda)
y <- c(1, 3, 2, 12, 1, 1)
x <- c(33, 14, 27, 90, 12, 17)


# a --------------------------------------------------
# Y_i represents the observed disease rate of i th county
# theta_i represents the expected disease rate of i th county


# b --------------------------------------------------
# Thetas are independent of each other.
# Thus, p(theta[i] | theta[-i], a, b, x, y) = p(theta[i] | a, b, xi, yi)

# sample_theta <- \(a, b, xi, yi) rgamma(1, yi + a, xi + b)

# all theta can be sampled in this problem setting (on R)
# x, y contain all observations (6 values here)
sample_theta <- \(a, b, x, y) rgamma(length(y), y + a, x + b)


# c --------------------------------------------------

# ap, bp are proposed values.
# as, bs are the current states of a, b.
# lr_ab <- \(ap, bp, as, bs, theta) {
#   return(
#     6 * (ap * log(bp) - as * log(bs) + lgamma(as) - lgamma(ap)) +
#       sum((ap - as) * log(theta) - (bp - bs) * theta) +
#       9 * (log(bp) - log(bs)) -
#       (ap - as + bp - bs)
#   )
# }

# But, the concise function below works
lr_ab <- \(ap, bp, as, bs, theta) {
  return(
    sum(dgamma(theta, ap, bp, log = TRUE) - dgamma(theta, as, bs, log = TRUE)) +
      dgamma(ap, 1, 1, log = TRUE) - dgamma(as, 1, 1, log = TRUE) +
      dgamma(bp, 10, 1, log = TRUE) - dgamma(bs, 10, 1, log = TRUE)
  )
}


# d --------------------------------------------------

sample_ab <- \(as, bs, theta, d1, d2) {
  ap <- abs(runif(1, as - d1, as + d1))
  bp <- abs(runif(1, bs - d2, bs + d2))

  if (log(runif(1)) < lr_ab(ap, bp, as, bs, theta)) {
    return(list(a = ap, b = bp, acs = 1))
  }
  return(list(a = as, b = bs, acs = 0))
}

mhmcmc <- \(y, x, d1, d2, n_iter) {
  post_samples <- list(
    THETA = matrix(nrow = n_iter, ncol = length(y)),
    A = numeric(n_iter),
    B = numeric(n_iter),
    ACS = 0
  )

  # init
  a <- 1
  b <- 10
  theta <- y / x

  for (i in seq_len(n_iter)) {
    theta <- sample_theta(a, b, x, y)
    ab_acs <- sample_ab(a, b, theta, d1, d2)
    a <- ab_acs$a
    b <- ab_acs$b

    post_samples$THETA[i, ] <- theta
    post_samples$A[i] <- a
    post_samples$B[i] <- b
    post_samples$ACS <- post_samples$ACS + ab_acs$acs
  }
  post_samples$ACS <- post_samples$ACS / n_iter
  return(post_samples)
}

n_iter <- 20000
post_samples <- mhmcmc(y, x, 1.5, 8, n_iter)
post_samples$ACS

samples_thin <- (\(mcmc_sample) {
  lag <- 10
  thinning <- seq(lag, n_iter, by = lag)
  return(list(
    THETA = mcmc_sample$THETA[thinning, ],
    A = mcmc_sample$A[thinning],
    B = mcmc_sample$B[thinning]
  ))
})(post_samples)

str(samples_thin)

acf(samples_thin$A)
acf(samples_thin$B)

coda::effectiveSize(samples_thin$A)
coda::effectiveSize(samples_thin$B)

par(mfrow = c(2, 1))
plot(samples_thin$A, type = "l")
plot(samples_thin$B, type = "l")
dev.off()


# e --------------------------------------------------

# i
png("fig/11-5-e-i.png", width = 900, height = 600)
par(mfrow = c(2, 3))
for (i in 1:6) {
  plot(density(samples_thin$THETA[, i], adjust = 2),
    main = bquote(theta[.(i)]),
    xlim = range(samples_thin$THETA)
  )
  abline(v = y[i] / x[i], lty = "dashed")
}
dev.off()

# ii
ab_prior <- rgamma(2000, 1, 1) / rgamma(2000, 10, 1)
png("fig/11-5-e-ii.png", width = 600, height = 600, res = 100)
plot(with(samples_thin, density(A / B, adjust = 2)), main = "a / b", xlab = "")
lines(density(ab_prior, adjust = 2), lty = "dashed")
abline(v = mean(y / x), lwd = 2)
legend("topright",
  legend = c("Prior", "Posterior", "mean(y / x)"),
  lty = c("dashed", "solid", "solid"), lwd = c(1, 1, 2)
)
dev.off()

# iii
theta2 <- samples_thin$THETA[, 2]
range_theta <- range(samples_thin$THETA)
png("fig/11-5-e-iii.png", width = 900, height = 600, res = 100)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2, 1))
for (i in 1:6) {
  plot(theta2, samples_thin$THETA[, i],
    xlim = range_theta, ylim = range_theta,
    xlab = expression(theta[2]), ylab = bquote(theta[.(i)])
  )
  abline(0, 1)
}
dev.off()

js <- c(1, 3, 4, 5, 6)
prob2overj <- apply(samples_thin$THETA[, js], 2, \(x) mean(theta2 > x))
names(prob2overj) <- c(1, 3, 4, 5, 6)
prob2overj

prob2max <- with(samples_thin, table(apply(THETA, 1, which.max)) / nrow(THETA))
prob2max

y / x
