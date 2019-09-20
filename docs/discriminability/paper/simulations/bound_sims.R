##----------------------------------------------
# Bound Simulation
##----------------------------------------------

require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(lolR)
require(abind)
require(emdbook)
no_cores = detectCores() - 1


## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  ni <- rowSums(rmultinom(n, 1, prob=rep(1/K, K)))
  X <- do.call(rbind, lapply(1:K, function(k) mvrnorm(n=ni[k], mus[,k], Sigmas[,,k])))
  Y <- do.call(c, sapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X, Y=Y))
}

sim.no_signal.bayes <- function(n=10000, d, sigma=1, npts=100) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  samp$X <- samp$X + array(rnorm(n*d), dim=c(n, d))*sigma
  samp$Z <- samp$Y

  range.x <- c(min(samp$X[,1]), max(samp$X[,1]))
  range.y <- c(min(samp$X[,2]), max(samp$X[,2]))

  probs <- sapply(1:2, function(y) mean(samp$Z == y))
  test.first=kde2d(samp$X[samp$Z == 1,1], samp$X[samp$Z == 1,2], lims=c(range.x, range.y), n=npts)
  test.first$z <- test.first$z/sum(test.first$z)
  test.sec=kde2d(samp$X[samp$Z == 2,1], samp$X[samp$Z == 2,2], lims=c(range.x, range.y), n=npts)
  test.sec$z <- test.sec$z/sum(test.sec$z)

  integration_reg <- array(NaN, dim=c(npts, npts))

  bayes <- sum(pmin(test.first$z, test.sec$z))*0.5

  return(bayes)
}

## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  ni <- rowSums(rmultinom(n, 1, prob=rep(1/K, K)))
  X <- do.call(rbind, lapply(1:K, function(k) mvrnorm(n=ni[k], mus[,k], Sigmas[,,k])))
  Y <- do.call(c, sapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X, Y=Y))
}

## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n, d, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n, priors=c(0.5,0.5))
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y, Z=samp$Y))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n, d, n.bayes=10000, n.pts=100,sigma=0) {
  S.class <- diag(d)
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus.class <- mvrnorm(n=2, c(0,0), S.class)
  pi.k <- 0.5  # equal chance of a new sample being from class 1 or class 2
  samp <- mgc:::mgc.sims.sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n,
                                 priors=c(pi.k, pi.k))
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y, Z=samp$Y))
}

## Crossed Signal Difference
# a simulation where classes are crossed but distinguishable
# also contains correlation btwn dimensions
# 2 classes
sim.crossed_sig <- function(n, d, K=16, sigma=0) {
  # class mus
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus.class <- t(mvrnorm(n=K, mu.class, S.class))
  ni <- n/K

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))  # covariances are orthogonal
  mus=cbind(rep(0, d), rep(0, d))

  X <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    rho <- runif(1, min=-1, max=1)*sqrt(sum(diag(Sigma.1)))
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rho
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rho
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- mgc:::mgc.sims.sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas,
                                  ni, priors=c(0.5, 0.5))
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X=X, Y=Y))
}

## Samples from Multiclass Gaussians
# a simulation where there are multiple classes present, and a correlation structure
# 2 classes
sim.multiclass_gaussian <- function(n, d, K=16, sigma=0) {
  # centers for the individuals
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  # sample the individual centers
  mus.class <- t(mvrnorm(n=K, mu.class, S.class))

  # individual covariances
  S.k <- diag(d)
  S.k[upper.tri(S.k)] <- 0.5  # correlated
  S.k[lower.tri(S.k)] <- 0.5
  Sigmas <- abind(lapply(1:K, function(k) S.k), along=3)

  # sample individuals w.p. 1/K
  samp <- mgc:::mgc.sims.sim_gmm(mus=mus, Sigmas=Sigmas, n, priors=rep(1/K, K))
  return(list(X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=samp$Y))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc <- function(n, d, K=16, sigma=0) {
  # centers
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class, S.class)))

  # probability of being each individual is 1/K
  ni <- rowSums(rmultinom(n, 1, prob=rep(1/K, K)))

  # individuals are either (1) a ball, or (2) a disc, around means
  X <- do.call(rbind, lapply(1:(K/2), function(k) {
    n.ball <- ni[2*(k-1)+1]; n.disc <- ni[2*k]
    X <- array(NaN, dim=c((n.ball + n.disc), d))
    X[1:n.ball,] <- sweep(mgc.sims.2ball(n.ball, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(n.ball + 1):(n.ball + n.disc),] <- sweep(mgc.sims.2sphere(n.disc, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=Y))
}

