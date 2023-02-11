if (!require(coda)) {
  install.packages("coda")
  library(coda)
}
dat <- scan("data/glucose.dat")


# a --------------------------------------------------

hist(dat, freq = FALSE)
lines(density(dat))

# b --------------------------------------------------
# see the djoint_gibbs function in section c

# c --------------------------------------------------

djoint_gibbs <- \(y, init_X, init_p, init_the_1, init_the_2,
  init_prec_1, init_prec_2, a, b, mu0, tau2_0, nu0, sig2_0, n_iter = 10000) {
  ny <- length(y)

  Xs <- matrix(nrow = n_iter, ncol = ny)
  Xs[1, ] <- X <- init_X
  PHI <- matrix(nrow = n_iter, ncol = 5)
  PHI[1, ] <- phi <- c(
    init_p, init_the_1, init_the_2, init_prec_1, init_prec_2
  )
  colnms <- c("p", "theta_1", "theta_2", "prec_1", "prec_2")
  colnames(PHI) <- names(phi) <- colnms

  vdnorm <- \(x, mu, sd) vapply(x, dnorm, numeric(1), mean = mu, sd = sd)
  for (s in 2:n_iter) {
    # x1, ..., xn
    props_1 <- phi["p"] * vdnorm(y, phi["theta_1"], sqrt(1 / phi["prec_1"]))
    props_2 <- (1 - phi["p"]) * vdnorm(y, phi["theta_2"], sqrt(1 / phi["prec_2"]))
    probs_1 <- props_1 / (props_1 + props_2)
    X <- sapply(probs_1, \(x) sample(c(1, 2), 1, prob = c(x, 1 - x)))

    # p
    nx1 <- length(X[X == 1])
    phi["p"] <- rbeta(1, a + nx1, b + ny - nx1)

    # theta_1
    y1 <- y[X == 1]
    ny1 <- length(y1)
    tau2_1_post <- 1 / (1 / tau2_0 + ny1 * phi["prec_1"])
    mu1_post <- (mu0 / tau2_0 + sum(y1) * phi["prec_1"]) * tau2_1_post
    phi["theta_1"] <- rnorm(1, mu1_post, sqrt(tau2_1_post))

    # prec_1
    my1 <- mean(y1)
    nu1_post <- nu0 + ny1
    sig2_1_post <- (nu0 * sig2_0 + sum((y1 - my1)^2) +
      ny1 * (my1 - phi["theta_1"])^2) / nu1_post
    phi["prec_1"] <- rgamma(1, nu1_post / 2, nu1_post * sig2_1_post / 2)
    # print(c(nu1_post, sig2_1_post))

    # theta_2
    y2 <- y[X == 2]
    ny2 <- length(y2)
    tau2_2_post <- 1 / (1 / tau2_0 + ny2 * phi["prec_2"])
    mu2_post <- (mu0 / tau2_0 + sum(y2) * phi["prec_2"]) * tau2_2_post
    phi["theta_2"] <- rnorm(1, mu2_post, sqrt(tau2_2_post))

    # prec_2
    my2 <- mean(y2)
    nu2_post <- nu0 + ny2
    sig2_2_post <- (nu0 * sig2_0 + sum((y2 - my2)^2) +
      ny2 * (my2 - phi["theta_2"])^2) / nu2_post
    phi["prec_2"] <- rgamma(1, nu2_post / 2, nu2_post * sig2_2_post / 2)
    # print(c(nu2_post, sig2_2_post))

    PHI[s, ] <- phi
    Xs[s, ] <- X
    # print(phi)
  }
  return(list(X = Xs, PHI = PHI))
}

a <- b <- 1
mu0 <- 120
tau2_0 <- 200
sig2_0 <- 1000
nu0 <- 10

init_X <- sample(c(1, 2), length(dat), replace = TRUE, prob = c(.5, .5))
init_p <- 0.5
init_the_1 <- mean(dat[init_X == 1])
init_the_2 <- mean(dat[init_X == 2])
init_prec_1 <- 1 / var(dat[init_X == 1])
init_prec_2 <- 1 / var(dat[init_X == 2])

gibbs_sample <- djoint_gibbs(
  dat, init_X, init_p, init_the_1, init_the_2, init_prec_1, init_prec_2,
  a, b, mu0, tau2_0, nu0, sig2_0
)

theta_max <- apply(gibbs_sample$PHI[, c("theta_1", "theta_2")], 1, max)
theta_min <- apply(gibbs_sample$PHI[, c("theta_1", "theta_2")], 1, min)

png("fig/6-2-c.png", width = 800, height = 400)
par(mfrow = c(1, 2))
acf(theta_max)
acf(theta_min)
dev.off()

coda::effectiveSize(theta_max)
coda::effectiveSize(theta_min)


# d --------------------------------------------------

get_pred <- \(x, phi) {
  pred <- NULL
  if (x == 1) pred <- rnorm(1, phi["theta_1"], sqrt(1 / phi["prec_1"]))
  if (x == 2) pred <- rnorm(1, phi["theta_2"], sqrt(1 / phi["prec_2"]))
  return(pred)
}

PRED <- matrix(nrow = nrow(gibbs_sample$X), ncol = ncol(gibbs_sample$X))
for (i in seq_len(nrow(gibbs_sample$X))) {
  X <- gibbs_sample$X[i, ]
  phi <- gibbs_sample$PHI[i, ]
  PRED[i, ] <- sapply(X, \(x) get_pred(x, phi))
}

png("fig/6-2-d.png", res = 100)
par(mar = c(5, 4.5, 2, 2))
plot(density(dat, adjust = 1), ylim = c(0, 0.015), xlab = "glucose", main = "")
lines(density(PRED, adjust = 1), lty = "dashed")
legend("topright", legend = c("Data", "Prediction"), lty = c("solid", "dashed"))
dev.off()

# plot(density(dat, adjust = 2), ylim = c(0, 0.015))
# lines(density(PRED, adjust = 2), lty = "dashed")
