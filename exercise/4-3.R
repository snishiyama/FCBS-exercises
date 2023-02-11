y_a <- c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6)
y_b <- c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7)

scaled_mean <- function(x) mean(x) / sd(x)

calc_t_pred <- function(n, theta) scaled_mean(rpois(n, theta))


# a --------------------------------------------------

t_pred_a <- sapply(rgamma(1000, 237, 20), \(x) calc_t_pred(10, x))
(t_obs_a <- scaled_mean(y_a))
hist(t_pred_a, main = NULL, xlab = expression(t(y[a])))
abline(v = t_obs_a, lwd = 3)
box()


# b --------------------------------------------------

t_pred_b <- sapply(rgamma(1000, 125, 14), \(x) calc_t_pred(13, x))
(t_obs_b <- scaled_mean(y_b))
hist(t_pred_b, main = NULL, xlab = expression(t(y[b])))
abline(v = t_obs_b, lwd = 3)
box()


# save figures in a single file
png("fig/4-3.png")
par(mfrow = c(2, 1), mar = c(4.5, 5, 2, 2))
hist(t_pred_a, main = NULL, xlab = expression(t(y[a])))
abline(v = t_obs_a, lwd = 3)
box()
hist(t_pred_b, main = NULL, xlab = expression(t(y[b])))
abline(v = t_obs_b, lwd = 3)
box()
dev.off()
