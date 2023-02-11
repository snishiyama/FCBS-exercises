bach <- scan("data/menchild30bach.dat")
nobach <- scan("data/menchild30nobach.dat")

table(c(bach, nobach))

# a --------------------------------------------------
# p(y | theta) ~ multinomial distribution(n, theta)
# n is the total number of observations.
# see https://en.wikipedia.org/wiki/Multinomial_distribution

# R has dmutinom() for PMF of multinomial dist.

# b --------------------------------------------------
# skip

# c --------------------------------------------------

rdir <- \(nsamp = 1, a) {
  Z <- matrix(rgamma(length(a) * nsamp, a, 1), nsamp, length(a), byrow = TRUE)
  return(Z / apply(Z, 1, sum))
}

sample_post <- \(dat, a_prior, n_sample) {
  classes <- seq_along(a_prior) - 1 # 0 should be included
  n_per_cls <- table(c(dat, classes)) - 1
  a_post <- n_per_cls + a_prior
  return(rdir(n_sample, a_post))
}

(n_cls <- max(c(bach, nobach)) + 1) # +1 for class `0`
n_sample <- 10000
res_bach <- sample_post(bach, rep(1, n_cls), n_sample)
res_nobach <- sample_post(nobach, rep(1, n_cls), n_sample)

apply(res_bach, 2, mean)
apply(res_nobach, 2, mean)

# prediction

(E_YApred <- apply(res_bach, 2, mean) * length(bach))
(E_YBpred <- apply(res_nobach, 2, mean) * length(nobach))

sample_pred <- \(n_obs, prob) c(rmultinom(1, size = n_obs, prob = prob))

bach_pred <- apply(res_bach, 1, \(p) sample_pred(length(bach), p))
nobach_pred <- apply(res_nobach, 1, \(p) sample_pred(length(nobach), p))

cls <- seq_len(n_cls) - 1

table(c(nobach, cls)) - 1

png("fig/12-4-c.png", width = 800, height = 600, res = 100)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1))
plot(cls + 0.1, apply(bach_pred, 1, mean),
  type = "h", lwd = 2,
  xlim = c(-0.5, 6.5), ylim = c(0, 30),
  xlab = "# of children", ylab = "Frequency"
)
lines(cls - 0.1, table(c(bach, cls)) - 1, type = "h", lwd = 2, col = "blue")
legend("topright",
  legend = c(bquote(Y[A]), bquote(tilde(Y)[A])),
  lty = "solid", col = c("blue", "black")
)
plot(cls + 0.1, apply(nobach_pred, 1, mean),
  type = "h", lwd = 2,
  xlim = c(-0.5, 6.5), ylim = c(0, 80),
  xlab = "# of children", ylab = "Frequency"
)
lines(cls - 0.1, table(c(nobach, cls)) - 1, type = "h", lwd = 2, col = "blue")
legend("topright",
  legend = c(bquote(Y[B]), bquote(tilde(Y)[B])),
  lty = "solid", col = c("blue", "black")
)
dev.off()

# d --------------------------------------------------

cnt_bach <- table(c(bach, cls)) - 1
cnt_nobach <- table(nobach)

png("fig/12-4-d.png", width = 960, height = 480, res = 100)
par(mfrow = c(1, 2))
# index 1 represents zero children
plot(t(bach_pred[2:3, ]),
  col = "gray", bg = "gray", pch = 21,
  xlab = "One child", ylab = "Two children", main = "bach"
)
points(cnt_bach[2], cnt_bach[3], col = "red", bg = "red", pch = 21)
plot(t(nobach_pred[2:3, ]),
  col = "gray", bg = "gray", pch = 21,
  xlab = "One child", ylab = "Two children", main = "nobach"
)
points(cnt_nobach[2], cnt_nobach[3], col = "red", bg = "red", pch = 21)
dev.off()
