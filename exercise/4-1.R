# posterior dist for second survey is Beta(31, 21)

# random sampling
set.seed(1)
rtheta1 <- rbeta(5000, 58, 44)
rtheta2 <- rbeta(5000, 31, 21)

# Pr(theta_1 < theta_2 | prior, dat)
print(sprintf(
      "Pr(theta_1 < theta_2 | prior, data) = %.3f",
      mean(rtheta1 < rtheta2)
))


# visualize to check the result above
par(mfrow = c(2, 1))
# mathematical
par(mar = c(5, 5, 2, 2))
curve(dbeta(x, 58, 44), from = 0, to = 1)
curve(dbeta(x, 31, 21), from = 0, to = 1, add = TRUE)
# mc results
plot(density(rtheta1, adjust = 2), xlim = c(0, 1))
lines(density(rtheta2, adjust = 2))
dev.off()
