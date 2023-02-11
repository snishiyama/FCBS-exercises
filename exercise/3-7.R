n1 <- 15
a1 <- 2


# a --------------------

poster_a <- function(theta) dbeta(theta, 3, 14)
e_post_a <- 3 / (3 + 14)
mode_post_a <- 2 / (2 + 13)
var_post_a <- 3 / 17 * (14 / 17) * (1 / 18)

print(sprintf(
     "E = %.3f; Mode = %.3f; SD = %.3f",
     e_post_a, mode_post_a, sqrt(var_post_a)
))

png("fig/3-7-a.png", width = 600, height = 600, res = 100)
curve(poster_a(x),
     from = 0, to = 1,
     xlab = expression(theta),
     ylab = expression(paste(p, "(", theta, "|", y, ")"))
)
dev.off()

# c --------------------
pred_b <- function(y2) {
     # using exp and lgamma to avoid overflow when computing gamma()
     gamma_1 <- exp(lgamma(17) - (lgamma(3) + lgamma(14)))
     gamma_2 <- exp(lgamma(y2 + 3) + lgamma(292 - y2) - lgamma(295))
     return(choose(278, y2) * gamma_1 * gamma_2)
}

y2_range <- 0:278
png("fig/3-7-c.png", width = 600, height = 600, res = 100)
plot(y2_range, pred_b(y2_range),
     type = "h",
     xlab = expression(y[2]),
     ylab = expression(Pr(Y[2] == y[2] ~ "|" ~ Y[1] == 2))
)
dev.off()

e_y2 <- sum(y2_range * pred_b(y2_range))
var_y2 <- sum(y2_range^2 * pred_b(y2_range)) - e_y2^2
sd_y2 <- sqrt(var_y2)

print(sprintf("E = %.3f; SD = %.3f", e_y2, sd_y2))


# d --------------------
likeli_d <- function(y2) dbinom(y2, 278, 2 / 15)

png("fig/3-7-d.png", width = 600, height = 600, res = 100)
par(mar = c(5, 5, 2, 2))
plot(y2_range, likeli_d(y2_range),
     type = "h",
     xlab = expression(y[2]),
     ylab = expression(Pr(Y[2] == y[2] ~ "|" ~ hat(theta) == 2 / 15))
)
dev.off()
