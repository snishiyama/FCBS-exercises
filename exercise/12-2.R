dat <- read.table("data/azdiabetes_alldata.dat", header = TRUE)
str(dat)


# sampling functions --------------------------------------------------

rtnorm <- \(n, mu, sigma, a, b) {
  bound <- \(x) pnorm(x, mu, sigma) # equal to pnorm((x- mu) / sigma, 0, 1)
  u <- runif(n, bound(a), bound(b))
  return(qnorm(u, mu, sigma))
}

rwish_ <- \(nu, S) {
  p <- ncol(S)
  Z <- matrix(rnorm(nu * p), nrow = nu, ncol = p)
  Y <- Z %*% chol(S)
  return(t(Y) %*% Y)
}

# Z is n x p matrix (n is # of obs. p is # of variables)
# nu0, S0 is prior for Sigma
sample_Sigma <- \(Z, nu0, S0) {
  n <- nrow(Z)
  nun <- nu0 + n
  iSn <- solve(S0 + t(Z) %*% Z)
  return(solve(rwish_(nun, iSn)))
  # return(solve(monomvn::rwish(nun, iSn)))
}

sample_Z <- \(Z, Sigma, Y) {
  n <- nrow(Z)
  p <- ncol(Z)
  Zp <- matrix(nrow = n, ncol = p)
  for (j in seq_len(p)) {
    Sjc <- Sigma[j, -j] %*% solve(Sigma[-j, -j])
    sz <- sqrt(Sigma[j, j] - Sjc %*% Sigma[-j, j])
    ez <- Z[, -j] %*% t(Sjc)
    for (i in seq_len(n)) {
      is_smaller <- Y[, j] < Y[i, j]
      is_larger <- Y[i, j] < Y[, j]
      a <- ifelse(sum(is_smaller, na.rm = TRUE) > 0, max(Z[is_smaller, j], na.rm = TRUE), -Inf)
      b <- ifelse(sum(is_larger, na.rm = TRUE) > 0, min(Z[is_larger, j], na.rm = TRUE), Inf)
      Zp[i, j] <- rtnorm(1, ez[i], sz, a, b)
    }
  }
  return(Zp)
}

gibbs <- \(Y, nu0, S0, n_iter, n_thin = max(1, floor(n_iter / 1000))) {
  n <- nrow(Y)
  p <- ncol(Y)
  n_samp <- floor(n_iter / n_thin)
  params <- list(
    SIGMA = array(dim = c(n_samp, p, p))
    # Z = array(dim = c(n_samp, n, p)) # too large data size to save
  )
  sample_Sigma_ <- \(Z) sample_Sigma(Z, nu0, S0)

  # init
  Yr <- apply(Y, 2, rank, ties.method = "random", na.last = "keep")
  is_missing <- is.na(Yr)
  Sigma <- S0
  if (any(is_missing)) {
    Z <- matrix(nrow = n, ncol = p)
    for (j in 1:p) {
      is_miss_j <- is_missing[, j]
      n_obs <- sum(!is_miss_j)
      Z[!is_miss_j, j] <- qnorm(Yr[!is_miss_j, j] / (n_obs + 1))
      Z[is_miss_j, j] <- mean(Z[!is_miss_j, j])
    }
  } else {
    Z <- apply(Yr, 2, \(y) qnorm(y / (n + 1)))
  }

  start <- Sys.time()
  for (i in seq_len(n_iter)) {
    if (i %% 500 == 0) {
      elapsed <- sprintf("%.2f", difftime(Sys.time(), start, unit = "secs"))
      cat("Running", i, "th loop. Elapsed", elapsed, "sec\n")
    }
    Sigma <- sample_Sigma_(Z)
    Z <- sample_Z(Z, Sigma, Yr)

    if (i %% n_thin == 0) {
      s <- floor(i / n_thin)
      params$SIGMA[s, , ] <- Sigma
      # params$Z[s, , ] <- Z # too large data size to save
    }
  }
  return(params)
}


# some functions for transformation --------------------------------------------

cov2cor <- \(Sigma) {
  sdiag <- matrix(0, nrow = nrow(Sigma), ncol = ncol(Sigma))
  diag(sdiag) <- 1 / sqrt(diag(Sigma))
  return(sdiag %*% Sigma %*% sdiag)
}

cor2beta <- \(PSI) {
  p <- ncol(PSI)
  BETA <- matrix(nrow = p, ncol = p)
  for (j in seq_len(p)) {
    BETA[j, -j] <- PSI[j, -j] %*% solve(PSI[-j, -j])
  }
  diag(BETA) <- 0
  return(BETA)
}

is_not_cross_zero <- \(Beta) {
  ci <- quantile(Beta, c(0.025, 0.975))
  if (all(ci > 0)) {
    return("+")
  } else if (all(ci < 0)) {
    return("-")
  }
  return("")
}

transf <- \(targ, fn) {
  out <- array(dim = dim(targ))
  for (i in seq_len(dim(targ)[1])) {
    out[i, , ] <- fn(targ[i, , ])
  }
  return(out)
}


# a --------------------------------------------------

n_iter <- 25000

# remove missing values
fnm_200 <- "data/12_2-a.rds"
if (file.exists(fnm_200)) {
  res_200 <- readRDS(fnm_200)
} else {
  p <- ncol(dat)
  # it took 2000 secs on Intel mac
  res_200 <- gibbs(dat[1:200, ], p + 2, diag(1, p, p), n_iter = n_iter)
  saveRDS(res_200, fnm_200)
}

str(res_200)

res_200$PSI <- transf(res_rmna$SIGMA, cov2cor)
res_200$BETA <- transf(res_rmna$PSI, cor2beta)
ePSI_rmna <- apply(res_200$PSI, c(2, 3), mean)
colnames(ePSI_rmna) <- rownames(ePSI_rmna) <- colnames(dat)
round(ePSI_rmna, 2)

# dependency graph
res_200$DEPS <- apply(res_rmna$BETA, c(2, 3), is_not_cross_zero)
colnames(res_200$DEPS) <- rownames(res_rmna$DEPS) <- colnames(dat)
res_200$DEPS

# b --------------------------------------------------

# including missing values
fnm_300 <- "data/12-2-b.rds"
if (file.exists(fnm_300)) {
  res_300 <- readRDS(fnm_300)
} else {
  p <- ncol(dat)
  # it took 1600 secs on Intel mac
  res_300 <- gibbs(dat[1:300, ], p + 2, diag(1, p, p), n_iter)
  saveRDS(res_300, fnm_300)
}

res_300$PSI <- transf(res_300$SIGMA, cov2cor)
res_300$BETA <- transf(res_300$PSI, cor2beta)
res_300$ePSI <- apply(res_300$PSI, c(2, 3), mean)
round(res_300$ePSI, 2)

# dependency graph
res_300$DEPS <- apply(res_300$BETA, c(2, 3), is_not_cross_zero)
colnames(res_300$DEPS) <- rownames(res_300$DEPS) <- colnames(dat)
res_300$DEPS


# plot results --------------------------------------------------

plot_psi_beta <- \(res, varnms) {
  np <- length(varnms)
  ori <- par(mfrow = c(np, np), mar = c(3, 3, 3, 1))
  for (i in seq_len(np)) {
    for (j in seq_len(np)) {
      if (i == j) {
        plot(NULL, xaxt = "n", yaxt = "n", bty = "n", xlim = c(-1, 1), ylim = c(-1, 1))
        text(0, 0, adj = 0.5, labels = varnms[i], cex = 1.5)
      } else if (j > i) {
        nm <- bquote(Psi[.(i) * "," * .(j)])
        plot(density(res$PSI[, i, j], adj = 2), main = nm, bty = "n")
        abline(v = mean(res$PSI[, i, j]), lwd = 2)
        ci <- quantile(res$PSI[, i, j], c(0.025, 0.975))
        abline(v = ci, lty = "dashed")
        box(col = ifelse(all(ci > 0) || all(ci < 0), "red", "black"))
      } else {
        nm <- bquote(beta[.(i) * "," * .(j)])
        plot(density(res$BETA[, i, j], adj = 2), main = nm, bty = "n")
        abline(v = 0, lty = "dashed", col = "blue")
        box(col = ifelse(res$DEPS[i, j] == "", "black", "red"))
      }
    }
  }
  par(ori)
}

png("fig/12-2-a.png", width = 1000, height = 1000, res = 100)
plot_psi_beta(res_200, colnames(dat))
dev.off()

png("fig/12-2-b.png", width = 1000, height = 1000, res = 100)
plot_psi_beta(res_300, colnames(dat))
dev.off()
