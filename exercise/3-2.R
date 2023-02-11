# Calculate a, b from theta and n
func_a <- function(theta, n) theta * n
func_b <- function(theta, n) (1 - theta) * n

# Calculate Pr(theta > 0.5 | sum(Y) = 57)
post_over_50 <- function(n, theta) {
  a <- func_a(theta, n)
  b <- func_b(theta, n)
  1 - pbeta(0.5, a + 57, b + 43)
}

theta <- seq(0.1, 0.9, by = 0.1)
n <- c(1, 2, 8, 16, 32)
probs <- outer(n, theta, post_over_50)

png("fig/3-2.png", width = 600, height = 600, res = 100)
contour(n, theta, probs,
  xlab = expression(n[0]),
  ylab = expression(theta[0]),
  nlevels = 20,
  labcex = 1.2
)
# use below if you like colored one
# filled.contour(n, theta, probs,
#         xlab = expression(n[0]),
#         ylab = expression(theta[0]),
#         nlevels = 20)
dev.off()
