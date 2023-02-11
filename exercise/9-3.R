dat <- read.table("data/crime.dat", header = TRUE)
head(dat)

X <- as.matrix(dat[, 2:ncol(dat)])
X <- cbind(matrix(1, nrow = nrow(X)), X)
colnames(X)[1] <- "intercept"
y <- as.matrix(dat[, "y", drop = FALSE])


# a --------------------------------------------------

# OLS

ols <- \(X, y) solve(t(X) %*% X) %*% t(X) %*% y
B_ols <- ols(X, y)


# Bayesian

# hyper parameter
mc <- \(y, X) {
  g <- length(y)
  nu0 <- 2
  s20 <- 1

  S <- 1000

  n <- nrow(X)
  p <- ncol(X)
  Hg <- (g / (g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
  SSRg <- t(y) %*% (diag(1, nrow = n) - Hg) %*% y
  # sigma
  s2 <- 1 / rgamma(S, (nu0 + n) / 2, (nu0 * s20 + SSRg) / 2)
  # beta
  Vb <- g / (g + 1) * solve(t(X) %*% X)
  Eb <- Vb %*% t(X) %*% y
  E <- matrix(rnorm(S * p, 0, sqrt(s2)), S, p)
  beta <- t(t(E %*% chol(Vb)) + c(Eb))

  return(list(beta = beta, sigma2 = s2))
}

beta <- mc(y, X)

mean_ci <- \(x) {
  m <- mean(x)
  names(m) <- "mean"
  ci <- quantile(x, prob = c(0.025, 0.975))
  return(c(m, ci))
}

apply(beta$beta, 2, mean_ci)
B_ols

cv_ols <- \(X_train, y_train, X_test, y_test) {
  beta_train <- ols(X_train, y_train)
  pred <- X_test %*% beta_train
  err <- mean((y_test - pred)^2)
  return(err)
}

cv_mc <- \(X_train, y_train, X_test, y_test) {
  post <- mc(y_train, X_train)
  Eb <- apply(post$beta, 2, mean)
  pred <- X_test %*% Eb
  err <- mean((y_test - pred)^2)
  return(err)
}

separate <- \(X, y) {
  n_dat <- nrow(X)
  ids <- seq_len(n_dat)
  train_id <- sample(ids, ceiling(n_dat / 2))
  test_id <- ids[-train_id]

  return(list(
    X_train = X[train_id, ],
    y_train = y[train_id, , drop = FALSE],
    X_test = X[test_id, ],
    y_test = y[test_id, , drop = FALSE]
  ))
}


# b, c --------------------------------------------------

errs <- matrix(nrow = 1000, ncol = 2)
colnames(errs) <- c("ols", "bayes")
for (i in 1:1000) {
  a <- separate(X, y)
  errs[i, 1] <- cv_ols(a$X_train, a$y_train, a$X_test, a$y_test)
  errs[i, 2] <- cv_mc(a$X_train, a$y_train, a$X_test, a$y_test)
}

apply(errs, 2, mean)
