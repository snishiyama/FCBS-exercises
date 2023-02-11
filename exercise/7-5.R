if (!require(monomvn)) {
  install.packages("monomvn")
  library(monomvn)
}
library(mvtnorm)

dat <- read.table("data/interexp.dat", header = TRUE)


# a --------------------------------------------------

theta_hat_A <- mean(dat$yA, na.rm = TRUE)
sig2_hat_A <- var(dat$yA, na.rm = TRUE)
theta_hat_B <- mean(dat$yB, na.rm = TRUE)
sig2_hat_B <- var(dat$yB, na.rm = TRUE)
rho_hat <- cor(dat, use = "complete.obs")[2, 1]

ans_a <- c(theta_hat_A, sig2_hat_A, theta_hat_B, sig2_hat_B, rho_hat)
names(ans_a) <- c("theta_hat_A", "sig2_hat_A", "theta_hat_B", "sig2_hat_B", "rho_hat")
print(ans_a)


# b --------------------------------------------------

dat_c <- dat
complete_B <- \(yA) {
  theta_hat_B + (yA - theta_hat_A) * rho_hat * sqrt(sig2_hat_B / sig2_hat_A)
}
complete_A <- \(yB) {
  theta_hat_A + (yB - theta_hat_B) * rho_hat * sqrt(sig2_hat_A / sig2_hat_B)
}
dat_c$yB <- with(dat_c, ifelse(is.na(yB), complete_B(yA), yB))
dat_c$yA <- with(dat_c, ifelse(is.na(yA), complete_A(yB), yA))

with(dat_c, t.test(yA, yB, paired = TRUE))


# c --------------------------------------------------

gibbs_sampling <- \(dat, n_iter = 1000) {
  n <- nrow(dat)
  p <- ncol(dat)
  # prior
  mu0 <- ybar0 <- apply(dat, 2, mean, na.rm = TRUE)
  nu0 <- p + 1
  dat_nonNA <- dat[apply(!is.na(dat), 1, all), ]
  S0 <- (t(dat_nonNA) - ybar0) %*% t(t(dat_nonNA) - ybar0) / n

  # init
  Phi <- solve(S0)
  O <- !is.na(dat)
  dat_full <- dat
  for (j in seq_len(ncol(dat))) {
    dat_full[is.na(dat_full[, j]), j] <- mean(dat_full[, j], na.rm = TRUE)
  }

  THETA <- SIGMA <- Y_MISS <- NULL
  # main loop
  for (s in 1:n_iter) {
    # theta
    # mun <- apply(dat_full, 2, mean) # identical to ybar
    Ln <- solve(solve(Phi) + n * Phi)
    ybar <- apply(dat_full, 2, mean)
    mun <- Ln %*% (solve(Phi) %*% mu0 + n * Phi %*% ybar)
    theta <- mvtnorm::rmvnorm(1, mun, Ln)

    # Sigma
    Sn <- S0 + (t(dat_full) - c(theta)) %*% t(t(dat_full) - c(theta))
    Phi <- monomvn::rwish(nu0 + n, solve(Sn))
    Sigma <- solve(Phi)

    # missing values
    for (i in 1:n) {
      a <- O[i, ] # non-NA cols
      b <- !O[i, ] # NA cols
      # skip rows where all variables are available
      if (all(a)) {
        next
      }
      iSa <- solve(Sigma[a, a]) # different from iSa <- Phi[a, a]
      beta_j <- Sigma[b, a] %*% iSa
      Sigma_j <- Sigma[b, b] - Sigma[b, a] %*% iSa %*% Sigma[a, b]
      theta_j <- theta[b] + beta_j %*% (t(dat_full[i, a]) - theta[a])
      dat_full[i, b] <- mvtnorm::rmvnorm(1, theta_j, Sigma_j)
    }
    THETA <- rbind(THETA, theta)
    SIGMA <- rbind(SIGMA, c(Sigma))
    Y_MISS <- rbind(Y_MISS, dat_full[!O])
  }
  return(list(THETA = THETA, SIGMA = SIGMA, Y_MISS = Y_MISS))
}

gibbs_sample <- gibbs_sampling(dat, 5000)

diff_theta_post <- apply(gibbs_sample$THETA, 1, \(x) x[1] - x[2])
plot(density(diff_theta_post, adjust = 2))
quantile(diff_theta_post, c(0.025, 0.50, 0.975))
mean(diff_theta_post)

acf(diff_theta_post) # auto-correlation
