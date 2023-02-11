rsize <- 10000


# a --------------------------------------------------

rtheta_a <- rgamma(rsize, 237, 20)
rtheta_b <- rgamma(rsize, 125, 14)

mean(rtheta_b < rtheta_a)


# b --------------------------------------------------

get_rtheta_b <- function(n0) rgamma(rsize, 12 * n0 + 113, n0 + 13)
ratio_a_over_b <- function(n0) mean(get_rtheta_b(n0) < rtheta_a)

ratios <- sapply(1:50, ratio_a_over_b)
plot(1:50, ratios,
     xlab = expression(n[0]),
     ylab = expression(Pr(theta[B] < theta[A] ~ "|" ~ y[A], y[B]))
)


# c --------------------------------------------------

pred_a <- rpois(rsize, rgamma(rsize, 237, 20))
pred_b <- rpois(rsize, rgamma(rsize, 125, 14))

mean(pred_b < pred_a)

get_pred_b <- function(n0) rpois(rsize, rgamma(rsize, 12 * n0 + 113, n0 + 13))
ratio_a_over_b_pred <- function(n0) mean(get_pred_b(n0) < pred_a)

ratios_pred <- sapply(1:50, ratio_a_over_b_pred)

par(mar = c(5, 5, 2, 2))
plot(1:50, ratios_pred,
     xlab = expression(n[0]),
     ylab = expression(Pr(tilde(Y)[B] < tilde(Y)[A] ~ "|" ~ y[A], y[B]))
)
dev.off()
