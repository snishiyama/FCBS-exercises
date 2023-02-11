dat_react <- read.table("data/cancer_react.dat", header = TRUE)
dat_noreact <- read.table("data/cancer_noreact.dat", header = TRUE)


# a --------------------------------------------------

# p(theta | data) is proportional to theta ^ (a + sum(yi) - 1) * exp(-theta * b + sum(xi))
# It follows gamma(a + sum(yi) , b + sum(xi))
# poster_prop <- function(theta, a, b, x, y) {
#   # theta ^ (a + sum(y) - 1) * exp(-theta * (b + sum(x)))
#   exp((sum(y) + a - 1) * log(theta) - theta * (b + sum(x)))
# }
poster <- function(theta, a, b, x, y) dgamma(theta, a + sum(y), b + sum(x))


# b --------------------------------------------------

poster_react <- function(x) poster(x, 1, 1, dat_react$x, dat_react$y)
poster_noreact <- function(x) poster(x, 1, 1, dat_noreact$x, dat_noreact$y)

curve(poster_react(x), from = 0, to = 10)
curve(poster_noreact(x), from = 0, to = 10)


# c ----------------------------------------

e_poster <- function(a, b, x, y) (a + sum(y)) / (b + sum(x))
ci_poster <- function(a, b, x, y) qgamma(c(0.025, 0.975), a + sum(y), b + sum(x))

for_opinion <- function(op_id, a1, a2, b1, b2, dat1, dat2) {
  e_post_1 <- e_poster(a1, b1, dat1$x, dat1$y)
  ci_post_1 <- ci_poster(a1, b1, dat1$x, dat1$y)
  e_post_2 <- e_poster(a2, b2, dat2$x, dat2$y)
  ci_post_2 <- ci_poster(a2, b2, dat2$x, dat2$y)
  print(sprintf(
    "E[theta1|dat1] = %.4f, 95CI = [%.4f, %.4f]",
    e_post_1, ci_post_1[1], ci_post_1[2]
  ))
  print(sprintf(
    "E[theta2|dat2] = %.4f, 95CI = [%.4f, %.4f]",
    e_post_2, ci_post_2[1], ci_post_2[2]
  ))

  mc_size <- 10000
  theta1_mc <- rgamma(mc_size, a1 + sum(dat1$y), b1 + sum(dat1$x))
  theta2_mc <- rgamma(mc_size, a2 + sum(dat2$y), b2 + sum(dat2$x))
  theta2_over_theta1 <- mean(theta2_mc > theta1_mc)
  print(sprintf("Pr(theta2 > theta1 | data) = %.4f", theta2_over_theta1))

  # plot
  dpost1 <- function(theta) poster(theta, a1, b1, dat1$x, dat1$y)
  dpost2 <- function(theta) poster(theta, a2, b2, dat2$x, dat2$y)
  thetas <- seq(0, 5, by = 0.01)
  plot(thetas, dpost1(thetas),
    type = "l", lwd = 2,
    main = sprintf("Opinion %d", op_id),
    xlab = expression(theta),
    ylab = expression(p(theta ~ "|" ~ data))
  )
  lines(thetas, dpost2(thetas), lwd = 2, col = "gray")
  legend(
    "topleft",
    legend = c(
      expression(p(theta[1] ~ "|" ~ data)),
      expression(p(theta[2] ~ "|" ~ data))
    ),
    bty = "n",
    lty = c("solid", "solid"),
    lwd = c(2, 2),
    col = c("black", "gray"),
    text.width = 1
  )
}

png("fig/4-5.png", width = 800, height = 800, res = 100)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
# opinion 1
for_opinion(1, 220, 220, 100, 100, dat_noreact, dat_react)
# opinion 2
for_opinion(2, 220, 2.2, 100, 1, dat_noreact, dat_react)
# opinion 3
for_opinion(3, 2.2, 2.2, 1, 1, dat_noreact, dat_react)
dev.off()
