dat_A <- scan("data/menchild30bach.dat")
dat_B <- scan("data/menchild30nobach.dat")

n_a <- length(dat_A)
n_b <- length(dat_B)

sum_a <- sum(dat_A)
sum_b <- sum(dat_B)

a_the_pre <- 2
b_the_pre <- 1


# b --------------------------------------------------

a_the_post <- sum_a + sum_b + a_the_pre
# b_the_post <- length(dat_A) + length(dat_B) * gamma + b_the_pre


# c --------------------------------------------------

# a_gam_post <- sum(dat_B) + a_gam_pre
# b_gam_post <- length(dat_B) * theta + b_gam_pre


# d --------------------------------------------------

djoint_gibbs <- function(gam_pre, init_t, init_g, n_iter = 5000) {
  a_gam_pre <- b_gam_pre <- gam_pre
  PHI <- matrix(nrow = n_iter, ncol = 2)
  PHI[1, ] <- phi <- c(init_t, init_g)

  a_gam_post <- sum_b + a_gam_pre
  for (s in 2:n_iter) {
    # theta
    b_the_post <- n_a + n_b * phi[2] + b_the_pre
    phi[1] <- rgamma(1, a_the_post, b_the_post)
    # gamma
    b_gam_post <- n_b * phi[1] + b_gam_pre
    phi[2] <- rgamma(1, a_gam_post, b_gam_post)
    PHI[s, ] <- phi
  }
  return(PHI)
}

gamma_pre <- c(8, 16, 32, 64, 128)
init_t <- mean(dat_A)
init_g <- mean(dat_B) / init_t

gibbs_sample <- lapply(gamma_pre, djoint_gibbs, init_t = init_t, init_g = init_g)
names(gibbs_sample) <- gamma_pre
e_diff <- sapply(gibbs_sample, \(phi) mean(phi[, 1] * (phi[, 2] - 1)))
e_diff
plot(gamma_pre, e_diff,
  type = "b",
  xlab = expression(a[gamma] == b[gamma]),
  ylab = expression(E(theta[B] - theta[A] ~ "|" ~ y[A], y[B]))
)
