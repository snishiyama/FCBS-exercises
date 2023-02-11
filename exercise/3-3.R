y_a <- c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6)
y_b <- c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7)


# a --------------------------------------------------

poster_a <- function(theta) dgamma(theta, 237, 20) # a = 120 + 117; b = 10 + 10
poster_b <- function(theta) dgamma(theta, 125, 14) # a = 12 + 113; b = 1 + 13
e_poster <- function(pri_a, pri_b, y) (pri_a + sum(y)) / (pri_b + length(y))
var_poster <- function(pri_a, pri_b, y) {
  (pri_a + sum(y)) / (pri_b + length(y))^2
}
ci_poster <- function(pri_a, pri_b, y, alpha = 0.05) {
  qgamma(c(alpha / 2, 1 - alpha / 2), pri_a + sum(y), pri_b + length(y))
}

ci_poster_a <- ci_poster(120, 10, y_a)
print(
  sprintf(
    "theta_a: E = %.4f, Var = %.4f, 95CI = [%.4f, %.4f]",
    e_poster(120, 10, y_a),
    var_poster(120, 10, y_a),
    ci_poster_a[1],
    ci_poster_a[2]
  )
)

ci_poster_b <- ci_poster(12, 1, y_b)
print(
  sprintf(
    "theta_a: E = %.4f, Var = %.4f, 95CI = [%.4f, %.4f]",
    e_poster(12, 1, y_b),
    var_poster(12, 1, y_b),
    ci_poster_b[1],
    ci_poster_b[2]
  )
)

png("fig/3-3-a.png")
curve(poster_a(x),
  from = 0, to = 20,
  xlab = expression(theta),
  ylab = expression(p(theta ~ "|" ~ y))
)
curve(poster_b(x), from = 0, to = 20, lty = "dashed", add = TRUE)
legend("topleft",
  legend = c(
    expression(p(theta[A] ~ "|" ~ y[A])),
    expression(p(theta[B] ~ "|" ~ y[B]))
  ),
  lty = c("solid", "dashed"),
  text.width = 3,
  y.intersp = 2
)
dev.off()


# b --------------------------------------------------

e_poster_b_n0 <- function(n0) e_poster(12 * n0, n0, y_b)

png("fig/3-3-b.png")
plot(n0 <- 1:50, e_poster_b_n0(n0),
  ylim = c(8.5, 12.5),
  xlab = expression(n[0]),
  ylab = expression(E * "[" * theta[B] ~ "|" ~ y[B] * "]")
)
abline(h = e_poster(120, 10, y_a), lty = "solid")
abline(h = ci_poster_a[1], lty = "dashed") # draw lower bound of ci
dev.off()


# c --------------------------------------------------

# b_post_fm_a_post <- function(theta) dgamma(theta, 350, 33)
# e_poster(237, 20, y_b)

# curve(b_post_fm_a_post(x),
#   from = 0, to = 20,
#   xlab = expression(theta),
#   ylab = expression(p(theta ~ "|" ~ y))
# )
# abline(v = e_poster(237, 20, y_b))
# curve(poster_a(x), from = 0, to = 20, lty = "dashed", add = TRUE)
# curve(poster_b(x), from = 0, to = 20, lty = "dotted", add = TRUE)
