n <- 43
y <- 15

e_poster <- function(a, b) (a + y) / (a + b + n) # posterior mean
mode_poster <- function(a, b) (a + y - 1) / (a + b + n - 2) # posterior mode
var_poster <- function(a, b) {
  numer <- e_poster(a, b) * (1 - e_poster(a, b))
  denom <- a + b + n + 1

  return(numer / denom)
}


# a --------------------------------------------------

prior_a <- function(theta) dbeta(theta, 2, 8)
likeli_a <- function(theta) dbinom(y, n, theta)
poster_a <- function(theta) dbeta(theta, 2 + y, 8 + n - y)

# ci of posterior distribution
poster_ci_a <- function(alpha) {
  lower <- alpha / 2
  upper <- 1 - alpha / 2
  qbeta(c(lower, upper), 2 + y, 8 + n - y)
}

png("fig/3-4-a.png", width = 1200, height = 1200, res = 100)
par(mfrow = c(2, 2))
curve(prior_a(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta))
)
curve(likeli_a(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(Pr(y ~ "|" ~ theta))
)
curve(poster_a(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta ~ "|" ~ y))
)
dev.off()

print(sprintf("E[theta | y] = %.4f", e_poster(2, 8)))
print(sprintf("Mode[theta | y] = %.4f", mode_poster(2, 8)))
print(sprintf("Var[theta | y] = %.4f", sqrt(var_poster(2, 8))))
poster_ci_a(0.05)


# b --------------------

prior_b <- function(theta) dbeta(theta, 8, 2)
poster_b <- function(theta) dbeta(theta, 8 + y, 2 + n - y)
poster_ci_b <- function(alpha) {
  qbeta(c(alpha / 2, 1 - alpha / 2), 8 + y, 2 + n - y)
}

png("fig/3-4-b.png", width = 1200, height = 800, res = 100)
par(mfrow = c(1, 2))
curve(prior_b(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta))
)
curve(poster_b(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta ~ "|" ~ y))
)
dev.off()

print(sprintf("E[theta | y] = %.4f", e_poster(8, 2)))
print(sprintf("Mode[theta | y] = %.4f", mode_poster(8, 2)))
print(sprintf("Var[theta | y] = %.4f", sqrt(var_poster(8, 2))))
poster_ci_b(0.05)


# c --------------------

beta_mix <- function(theta) {
  term_theta <- 3 * theta * (1 - theta)^7 + theta^7 * (1 - theta)
  term_scale <- gamma(10) / (gamma(2) * gamma(8))
  return(0.25 * term_scale * term_theta)
}

# This returns the same results
# beta_mix <- function(theta) {
#   return(0.75 * dbeta(theta, 2, 8) + 0.25 * dbeta(theta, 8, 2))
# }

png("fig/3-4-c.png", width = 600, height = 600, res = 100)
curve(beta_mix(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta))
)
dev.off()


# d iii --------------------

# p(theta) x p(y | theta)
poster_mix <- function(theta) {
  term_gamma <- gamma(10) / (gamma(2) * gamma(8))
  term_theta_1 <- theta^16 * (1 - theta)^35
  term_theta_2 <- theta^22 * (1 - theta)^29
  return(0.25 * term_gamma * choose(43, 15) * (3 * term_theta_1 + term_theta_2))
}

# png("fig/3-4-d.png", width = 600, height = 600, res = 100)
curve(poster_mix(x),
  from = 0, to = 1,
  xlab = expression(theta),
  ylab = expression(p(theta ~ "|" ~ y))
)
dev.off()

# find mode numerically
theta_range <- seq(0, 1, by = .0001)
theta_range[which.max(poster_mix(theta_range))]


# some extra code fragments --------------------------------------------------
# to figure out how to write the functions that return the same results
# of poster_mix() using dbeta()

# simply mix posteriors from the two priors dbeta(theta, 17, 36) and dbeta(theta, 23, 30)
# This returns different results
poster_mix2 <- function(theta) {
  0.75 * dbeta(theta, 17, 36) + 0.25 * dbeta(theta, 23, 30)
}

theta_range[which.max(poster_mix2(theta_range))]

# This pattern returns same results
poster_mix3 <- function(theta) {
  gamma_1 <- exp(lgamma(17) + lgamma(36))
  gamma_2 <- exp(lgamma(23) + lgamma(30))
  3 * gamma_1 * dbeta(theta, 17, 36) + gamma_2 * dbeta(theta, 23, 30)
}

theta_range[which.max(poster_mix3(theta_range))]

# posterior distribution (scaling the poster_mix3 distribution)
poster_mix_scaled <- function(theta) {
  numer <- poster_mix3(theta)
  gamma_1 <- exp(lgamma(17) + lgamma(36))
  gamma_2 <- exp(lgamma(23) + lgamma(30))
  denom <- 3 * gamma_1 + gamma_2
  return(numer / denom)
}

# Integral of poster_mix_scaled is 1
integrate(poster_mix_scaled, lower = 0, upper = 1)
