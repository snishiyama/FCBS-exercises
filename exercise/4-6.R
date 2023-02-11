log_odds <- function(theta) log(theta) - log(1 - theta)

set.seed(1)
rand_beta_prior <- rbeta(10000, 1, 1)
hist(log_odds(rand_beta_prior), freq = FALSE, main = NULL)
# lines(density(log_odds(rand_beta_prior)), add = TRUE)
box()
