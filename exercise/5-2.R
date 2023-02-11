calc_mun <- \(mu0, k0, y_bar, n) (k0 * mu0 + n * y_bar) / (k0 + n)

calc_sig2n <- function(mu0, k0, sig2_0, nu0, y_bar, s, n) {
  ss <- s^2 * (n - 1)
  numer <- nu0 * sig2_0 + ss + (k0 * n * (y_bar - mu0)^2) / (k0 + n)
  return(numer / (nu0 + n))
}

get_post_sample <- function(mu0, k0, s2_0, nu0, y_bar, s, n) {
  mc_size <- 10000
  s2n <- calc_sig2n(mu0, k0, s2_0, nu0, y_bar, s, n)
  nun <- nu0 + n
  s2n_post_mc <- 1 / rgamma(mc_size, nun / 2, s2n * nun / 2)

  mun <- calc_mun(mu0, k0, y_bar, n)
  kn <- k0 + n
  mu_post_mc <- rnorm(mc_size, mun, sqrt(s2n_post_mc / kn))

  return(mu_post_mc)
}

mu0 <- 75
sig2_0 <- 100

n_a <- n_b <- 16
y_bar_a <- 75.2
s_a <- 7.3
y_bar_b <- 77.5
s_b <- 8.1

get_post_sample_a <- \(k0, nu0) {
  get_post_sample(mu0, k0, sig2_0, nu0, y_bar_a, s_a, n_a)
}

get_post_sample_b <- \(k0, nu0) {
  get_post_sample(mu0, k0, sig2_0, nu0, y_bar_b, s_b, n_b)
}

vec_k0 <- vec_nu0 <- 1:32

post_sample_a <- mapply(get_post_sample_a, k0 = vec_k0, nu0 = vec_nu0)
post_sample_b <- mapply(get_post_sample_b, k0 = vec_k0, nu0 = vec_nu0)

probs <- apply(post_sample_a < post_sample_b, MARGIN = 2, FUN = mean)
names(probs) <- vec_k0
print(probs[c(1, 2, 4, 8, 16, 32)])
plot(vec_k0, probs)

png("fig/5-2.png", width = 600, height = 600, res = 100)
par(mfrow = c(2, 3), mar = c(3, 3, 4, 2))
for (i in c(1, 2, 4, 8, 16, 32)) {
  title <- bquote({
    kappa[0] == nu[0]
  } == .(as.character(i)) * ";"
  ~ Pr(theta[A] < theta[B]) == .(sprintf("%.3f", probs[i])))
  plot(specs <- density(post_sample_a[, i], adjust = 2), main = title)
  lines(density(post_sample_b[, i], bw = specs$bw), lty = "dashed")
}
dev.off()
