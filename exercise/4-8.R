bach <- scan("data/menchild30bach.dat")
nobach <- scan("data/menchild30nobach.dat")

mc_size <- 5000
post_mc <- \(a, b, dat) rgamma(mc_size, a + sum(dat), b + length(dat))

a_pre <- 1
b_pre <- 2

post_a <- post_mc(a_pre, b_pre, bach)
post_b <- post_mc(a_pre, b_pre, nobach)


# a --------------------------------------------------

pred_a <- rpois(mc_size, post_a)
pred_b <- rpois(mc_size, post_b)

tab_pred_a <- table(pred_a) / mc_size
tab_pred_b <- table(pred_b) / mc_size

png("fig/4-8-a.png", res = 100)
par(mar = c(5, 4.5, 2, 2))
plot(as.numeric(names(tab_pred_a)) - .1, as.vector(tab_pred_a),
     col = gray(0), type = "h", lwd = 4,
     xlab = "# of children", ylab = "Density"
)
points(as.numeric(names(tab_pred_b)) + .1, tab_pred_b,
     col = "gray", type = "h", lwd = 4
)
legend("topright",
     legend = c(expression(tilde(Y)[A]), expression(tilde(Y)[B])),
     lty = "solid", col = c("black", "gray")
)
dev.off()


# b --------------------------------------------------

quantile(post_b - post_a, c(0.025, 0.975))
quantile(pred_b - pred_a, c(0.025, 0.975))

hist(post_b - post_a)
hist(pred_b - pred_a)

plot(density(post_a, adjust = 2), xlim = c(0.4, 1.8), ylim = c(0, 5))
lines(density(post_b, adjust = 2), lty = "dashed")


# c --------------------------------------------------

tab_nobach <- table(nobach)
xvals <- as.numeric(names(tab_nobach))
plot(xvals - .1, as.vector(tab_nobach) / length(nobach),
     col = gray(0), type = "h", lwd = 4, xlim = c(0, 6), ylim = c(0, .4)
)
points(xvals + .1, dpois(xvals, 1.4),
     col = "gray", type = "h", lwd = 4
)


# d --------------------------------------------------

get_01 <- function(theta) {
     rsample <- rpois(218, theta)
     table(rsample)[c("0", "1")]
}

res <- vapply(post_b, get_01, numeric(2))

png("fig/4-8-d.png")
plot(res["0", ], res["1", ],
     col = clr <- rgb(0, 0, 0, alpha = 0.2), bg = clr,
     xlim = c(20, 90), ylim = c(40, 110), pch = 19,
     xlab = "# of men having zero children",
     ylab = "# of men having one child"
)
points(
     x = tab_nobach["0"], y = tab_nobach["1"],
     col = "red", bg = "red", pch = 88
)
dev.off()
