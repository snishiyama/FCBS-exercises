# a --------------------

# Beta dist with mode 1/3
prior_one_third <- function(theta, a) dbeta(theta, a, 2 * a - 1)
# Beta dist with mode 1/2
prior_half <- function(theta, a) dbeta(theta, a, a)

# a1, a2 are the parameter for Beta dists with mode 1/3 and 1/2, respectively
# ratio should be 0 to 1
prior_mix <- function(theta, a1, a2, ratio) {
  ratio * prior_one_third(theta, a1) + (1 - ratio) * prior_half(theta, a2)
}

# the equation below matches a + b between the two distributions
a1 <- 101
a2 <- (3 * a1 - 1) / 2

# 10 yen coin made in 1990
n1 <- 50
upside <- 26

# proportional to the posterior (i.e., the numerator of the bayes rule)
poster <- function(theta, a1, a2, ratio, n, k) {
  one_third <- (3 * a1 - 2) * choose(3 * a1 - 3, 2 * a1 - 2) *
    theta^(a1 + k - 1) * (1 - theta)^(2 * a2 - 1 + (n - k) - 1)
  half <- (2 * a2 - 1) * choose(2 * a2 - 2, a2 - 1) *
    theta^(a2 + k - 1) * (1 - theta)^(a2 + (n - k) - 1)
  return((ratio * one_third + (1 - ratio) * half) * choose(n, k))
}

png("fig/3-8.png", width = 1600, height = 800, res = 100)
par(mfrow = c(2, 4))
curve(prior_one_third(x, 3))
curve(prior_half(x, 4))
curve(prior_mix(x, 3, 4, 0.8))
curve(poster(x, 3, 4, 0.8, 50, 26))
curve(prior_one_third(x, 101))
curve(prior_half(x, 151))
curve(prior_mix(x, 101, 151, 0.8))
curve(poster(x, 101, 151, 0.8, 50, 26))
dev.off()
