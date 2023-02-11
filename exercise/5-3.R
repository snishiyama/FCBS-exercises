mu_n <- function(mu0, k0, dat) {
  n <- length(dat)
  return((k0 * mu0 + n * mean(dat)) / (k0 + n))
}

sig2_n <- function(mu0, k0, sig2_0, nu0, dat) {
  ss <- sum((dat - mean(dat))^2) # var(dat) * (n - 1)
  n <- length(dat)
  numer <- nu0 * sig2_0 + ss + (k0 * n * (mean(dat) - mu0)^2) / (k0 + n)
  return(numer / (nu0 + n))
}


get_post_sample <- function(dat) {
  mc_size <- 10000
  n <- length(dat)

  s2_n <- sig2_n(mu0, k0, sig2_0, nu0, dat)
  s2_post_mc <- 1 / rgamma(mc_size, (nu0 + n) / 2, s2_n * (nu0 + n) / 2)
  e_s2_n <- s2_n * (nu0 + n) * 0.5 / ((nu0 + n) * 0.5 - 1)

  mu_post <- mu_n(mu0, k0, dat)
  mu_post_mc <- rnorm(mc_size, mu_post, sqrt(s2_post_mc / (k0 + n)))

  return(list(
    mu_n = mu_post,
    e_s2_n = e_s2_n,
    mu_post_mc = mu_post_mc,
    s2_post_mc = s2_post_mc
  ))
}

# marginal posterior distribution of theta
mpost_theta <- function(theta, dat, mu0, sig2_0, nu0, k0) {
  n <- length(dat)
  kn <- k0 + n
  nu_n <- nu0 + n
  mun <- mu_n(mu0, k0, dat)
  sig2n <- sig2_n(mu0, k0, sig2_0, nu0, dat)

  term1 <- sqrt(kn / (2 * pi))
  term2_numer <- (nu_n * sig2n / 2)^(nu_n / 2)
  term2_denom <- exp(lgamma(nu_n / 2))
  term3_numer <- exp(lgamma((nu_n + 1) / 2))
  term3_denom <- (((kn * (theta - mun)^2) + nu_n * sig2n) / 2)^
    ((nu_n + 1) / 2)

  return(term1 * term2_numer / term2_denom * term3_numer / term3_denom)
}

# use 5-1 data and priors
dat_sc1 <- read.table("data/school1.dat", header = FALSE, col.names = "x")$x

mu0 <- 5
sig2_0 <- 4
k0 <- 1
nu0 <- 2

mpost_theta_sc1 <- \(x) mpost_theta(x, dat_sc1, mu0, sig2_0, nu0, k0)
post_sample_sc1 <- get_post_sample(dat_sc1)

x <- seq(6, 12, length.out = 100)
plot(x, mpost_theta_sc1(x), type = "l", col = "red")
lines(density(post_sample_sc1$mu_post_mc))
