if (!require(monomvn)) {
  install.packages("monomvn")
  library(monomvn)
}
library(mvtnorm)

blue <- read.table("data/bluecrab.dat", col.names = c("depth", "width"))
orange <- read.table("data/orangecrab.dat", col.names = c("depth", "width"))


# a --------------------------------------------------

gibbs_mvnorm <- function(dat, n_iter = 10000) {
  n <- nrow(dat)
  nu0 <- 4
  ybar <- apply(dat, 2, mean)
  S0 <- L0 <- Sigma <- cov(dat)
  THETA <- SIGMA <- NULL
  for (i in seq_len(n_iter)) {
    # theta
    Ln <- solve(solve(L0) + n * solve(Sigma))
    mun <- Ln %*% (solve(L0) %*% ybar + n * solve(Sigma) %*% ybar)
    theta <- mvtnorm::rmvnorm(1, mun, Ln)

    # sigma
    Y_minus_theta <- t(dat) - c(theta)
    Sn <- S0 + Y_minus_theta %*% t(Y_minus_theta)
    Sigma <- solve(monomvn::rwish(nu0 + n, solve(Sn)))

    THETA <- rbind(THETA, theta)
    SIGMA <- rbind(SIGMA, c(Sigma))
  }

  return(list(THETA = THETA, SIGMA = SIGMA))
}

param_gibbs_blue <- gibbs_mvnorm(blue)
param_gibbs_orange <- gibbs_mvnorm(orange)


# b --------------------------------------------------

png("fig/7-3-b.png")
par(mar = c(5, 4.5, 2, 2))
with(param_gibbs_blue, plot(THETA[, 1], THETA[, 2],
  col = "blue",
  xlim = c(10, 14), ylim = c(11, 17.5),
  xlab = expression(theta[1]), ylab = expression(theta[2])
))
with(param_gibbs_orange, points(THETA[, 1], THETA[, 2], col = "orange"))
dev.off()


# c --------------------------------------------------

head(param_gibbs_blue$SIGMA)

rho <- \(x) x[2] / (sqrt(x[1]) * sqrt(x[4]))

rho_blue <- apply(param_gibbs_blue$SIGMA, 1, rho)
rho_orange <- apply(param_gibbs_orange$SIGMA, 1, rho)

png("fig/7-3-c.png")
par(mar = c(5, 5, 2, 2))
plot(dens_blue <- density(rho_blue, adj = 2),
  xlim = c(0.9, 1.0), ylim = c(0, 120),
  main = "", xlab = expression(rho)
)
polygon(dens_blue, col = rgb(0, 0.5, 1, alpha = 0.5))
polygon(density(rho_orange, adj = 2), col = rgb(1, 0.5, 0, alpha = 0.5))
dev.off()

mean(rho_blue < rho_orange)

acf(rho_blue)
acf(rho_orange)


# extra --------------------------------------------------
# posterior predictive check

pred <- \(THETA, SIGMA) {
  Y <- matrix(nrow = nrow(THETA), ncol = 2)
  for (i in seq_len(nrow(THETA))) {
    Y[i, ] <- mvtnorm::rmvnorm(1, THETA[i, ], matrix(SIGMA[i, ], ncol = 2))
  }
  return(Y)
}

pred_blue <- with(param_gibbs_blue, pred(THETA, SIGMA))
pred_orange <- with(param_gibbs_orange, pred(THETA, SIGMA))


pars <- par(mfrow = c(1, 2))
plot(pred_blue[, 1], pred_blue[, 2], col = "gray", main = "Blue")
points(blue$depth, blue$width, pch = 88)

plot(pred_orange[, 1], pred_orange[, 2], col = "gray", main = "Orainge")
points(orange$depth, orange$width, pch = 88)
par(pars)
