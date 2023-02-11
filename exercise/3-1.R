# b --------------------

binary57 <- function(theta) {
  choose(100, 57) * theta^57 * (1 - theta)^43
}

theta_b <- seq(0, 1, by = 0.1)
val_b <- binary57(theta_b)

png("fig/3-1-b.png", width = 600, height = 600, res = 100)
plot(theta_b, val_b,
  type = "b",
  xlab = expression(theta),
  ylab = expression(Pr(Y == 57 ~ "|" ~ theta))
)
dev.off()


# c --------------------

numer_c <- val_b * (1 / length(theta_b))
denom_c <- sum(numer_c)
val_c <- numer_c / denom_c

png("fig/3-1-c.png", width = 600, height = 600, res = 100)
plot(theta_b, val_c,
  type = "b",
  xlab = expression(theta),
  ylab = expression(p(theta ~ "|" ~ Y == 57))
)
dev.off()


# d --------------------

theta_d <- seq(0, 1, by = 0.01)
val_d <- binary57(theta_d)

png("fig/3-1-d.png", width = 600, height = 600, res = 100)
plot(theta_d, val_d,
  type = "l",
  xlab = expression(theta),
  ylab = expression(Pr(Y == 57 ~ "|" ~ theta) ~ p(theta))
)
dev.off()


# e --------------------

png("fig/3-1-e.png", width = 600, height = 600, res = 100)
curve(dbeta(x, 58, 44), from = 0, to = 1)
dev.off()
