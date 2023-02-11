dat_sc1 <- read.table("data/school1.dat", header = FALSE, col.names = "time")
dat_sc2 <- read.table("data/school2.dat", header = FALSE, col.names = "time")
dat_sc3 <- read.table("data/school3.dat", header = FALSE, col.names = "time")

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

mu0 <- 5
sig2_0 <- 4
k0 <- 1
nu0 <- 2

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


# a --------------------------------------------------

print_ans_a <- function(school_id, post_sample) {
  mu_post <- post_sample$mu_n
  e_sd_n <- sqrt(post_sample$e_s2_n)
  mu_post_mc <- post_sample$mu_post_mc
  sd_post_mc <- sqrt(post_sample$s2_post_mc)

  ci_mu <- quantile(mu_post_mc, c(0.025, 0.975))
  ci_sd <- quantile(sd_post_mc, c(0.025, 0.975))

  ans_mu <- sprintf(
    "E[theta] = %.4f, 95CI = [%.4f, %.4f]",
    mu_post, ci_mu[1], ci_mu[2]
  )
  ans_sd <- sprintf(
    "E[sig] = %.4f, 95CI = [%.4f, %.4f]",
    e_sd_n, ci_sd[1], ci_sd[2]
  )
  print(school_id)
  print(ans_mu)
  print(ans_sd)
}

post_sample_sc1 <- get_post_sample(dat_sc1$time)
post_sample_sc2 <- get_post_sample(dat_sc2$time)
post_sample_sc3 <- get_post_sample(dat_sc3$time)

print_ans_a("School1", post_sample_sc1)
print_ans_a("School2", post_sample_sc2)
print_ans_a("School3", post_sample_sc3)


# b --------------------------------------------------

sc1_sc2 <- post_sample_sc1$mu_post_mc < post_sample_sc2$mu_post_mc
sc1_sc3 <- post_sample_sc1$mu_post_mc < post_sample_sc3$mu_post_mc
sc2_sc3 <- post_sample_sc2$mu_post_mc < post_sample_sc3$mu_post_mc
sc2_sc1 <- post_sample_sc2$mu_post_mc < post_sample_sc1$mu_post_mc
sc3_sc1 <- post_sample_sc3$mu_post_mc < post_sample_sc1$mu_post_mc
sc3_sc2 <- post_sample_sc3$mu_post_mc < post_sample_sc2$mu_post_mc

print(sprintf("1 < 2 < 3: %.4f", mean(sc1_sc2 & sc2_sc3)))
print(sprintf("1 < 3 < 2: %.4f", mean(sc1_sc3 & sc3_sc2)))
print(sprintf("2 < 1 < 3: %.4f", mean(sc2_sc1 & sc1_sc3)))
print(sprintf("2 < 3 < 1: %.4f", mean(sc2_sc3 & sc3_sc1)))
print(sprintf("3 < 1 < 2: %.4f", mean(sc3_sc1 & sc1_sc2)))
print(sprintf("3 < 2 < 1: %.4f", mean(sc3_sc2 & sc2_sc1)))


# c --------------------------------------------------

get_pred <- function(post_sample) {
  rnorm(10000, post_sample$mu_post_mc, sqrt(post_sample$s2_post_mc))
}

pred_sc1 <- get_pred(post_sample_sc1)
pred_sc2 <- get_pred(post_sample_sc2)
pred_sc3 <- get_pred(post_sample_sc3)

pred_sc1_sc2 <- pred_sc1 < pred_sc2
pred_sc1_sc3 <- pred_sc1 < pred_sc3
pred_sc2_sc1 <- pred_sc2 < pred_sc1
pred_sc2_sc3 <- pred_sc2 < pred_sc3
pred_sc3_sc1 <- pred_sc3 < pred_sc1
pred_sc3_sc2 <- pred_sc3 < pred_sc2

print(sprintf("1 < 2 < 3: %.4f", mean(pred_sc1_sc2 & pred_sc2_sc3)))
print(sprintf("1 < 3 < 2: %.4f", mean(pred_sc1_sc3 & pred_sc3_sc2)))
print(sprintf("2 < 1 < 3: %.4f", mean(pred_sc2_sc1 & pred_sc1_sc3)))
print(sprintf("2 < 3 < 1: %.4f", mean(pred_sc2_sc3 & pred_sc3_sc1)))
print(sprintf("3 < 1 < 2: %.4f", mean(pred_sc3_sc1 & pred_sc1_sc2)))
print(sprintf("3 < 2 < 1: %.4f", mean(pred_sc3_sc2 & pred_sc2_sc1)))


# d --------------------------------------------------

print(mean(sc2_sc1 & sc3_sc1))
print(mean(pred_sc2_sc1 & pred_sc3_sc1))

# check
# plots
par(mfrow = c(2, 1))
par(mar = c(5, 5, 2, 2))
dens_theta_sc1 <- density(post_sample_sc1$mu_post_mc, adj = 2)
plot(dens_theta_sc1,
  xlim = c(4, 12), lwd = 2, main = "",
  xlab = expression(theta)
)
lines(density(post_sample_sc2$mu_post_mc, bw = dens_theta_sc1$bw),
  lwd = 2, col = gray(.5)
)
lines(density(post_sample_sc3$mu_post_mc, bw = dens_theta_sc1$bw),
  lwd = 2, lty = "dashed"
)
legend("topleft",
  legend = c("School1", "School2", "School3"),
  lty = c("solid", "solid", "dashed"), lwd = c(2, 2, 2),
  col = c("black", gray(.5), "black"), text.width = 1.1, y.intersp = 2
)

dens_pred_sc1 <- density(pred_sc1, adj = 2)
plot(dens_pred_sc1, ylim = c(0, 0.15), lwd = 2, main = "", xlab = "Y")
lines(density(pred_sc2, bw = dens_pred_sc1$bw), lwd = 2, col = gray(.5))
lines(density(pred_sc3, bw = dens_pred_sc1$bw), lwd = 2, lty = "dashed")
legend("topleft",
  legend = c("School1", "School2", "School3"),
  lty = c("solid", "solid", "dashed"), lwd = c(2, 2, 2),
  col = c("black", gray(.5), "black"), text.width = 5, y.intersp = 2
)

dev.off()

# sum of probs is 1
sum(c(
  mean(sc1_sc2 & sc2_sc3),
  mean(sc1_sc3 & sc3_sc2),
  mean(sc2_sc1 & sc1_sc3),
  mean(sc2_sc3 & sc3_sc1),
  mean(sc3_sc1 & sc1_sc2),
  mean(sc3_sc2 & sc2_sc1)
))

sum(c(
  mean(pred_sc1_sc2 & pred_sc2_sc3),
  mean(pred_sc1_sc3 & pred_sc3_sc2),
  mean(pred_sc2_sc1 & pred_sc1_sc3),
  mean(pred_sc2_sc3 & pred_sc3_sc1),
  mean(pred_sc3_sc1 & pred_sc1_sc2),
  mean(pred_sc3_sc2 & pred_sc2_sc1)
))
