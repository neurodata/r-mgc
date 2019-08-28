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
no_cores = detectCores() - 1

## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n) {
  K <- dim(mus)[2]
  ni <- round(n/K)
  labs <- c(sapply(1:K, function(k) rep(k, ni)))
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs))
}


## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n, d, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d)), Y=samp$Y,
              Z=c(rep(1, n/2), rep(2, n/2))))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n, d, sigma=0) {
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus=cbind(rep(0, d), c(4, rep(0, d-1))) # with a mean signal shift between the classes
  samp <- sim_gmm(mus=mus, Sigmas=abind(Sigma, Sigma, along=3), n)
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y,
              Z=c(rep(1, n/2), rep(2, n/2))))
}

## Crossed Signal Difference
# a simulation where classes are crossed but distinguishable
# also contains correlation btwn dimensions
# 2 classes
sim.crossed_sig <- function(n, d, K=16, sigma=0) {
  # class mus
  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)*1.25
  S.class <- diag(d)*sqrt(K)

  mus.class <- t(rbind(mvrnorm(n=K/2, mu.class.1, S.class),
                       mvrnorm(n=K/2, mu.class.2, S.class)))
  ni <- n/K

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))
  mus=cbind(rep(0, d), rep(0, d))

  X <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    rho <- runif(1, min=-1, max=1)*sqrt(2*0.1)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rho
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rho
    sim <- sim_gmm(mus=mus, Sigmas=Sigmas, ni)
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X=X, Y=Y, Z=c(rep(1, n/2), rep(2, n/2))))
}

## Samples from Multiclass Gaussians
# a simulation where there are multiple classes present, and a correlation structure
# 2 classes
sim.multiclass_gaussian <- function(n, d, K=16, sigma=0) {
  S.k <- diag(d)*1
  S.k[upper.tri(S.k)] <- 0.5  # correlated
  S.k[lower.tri(S.k)] <- 0.5

  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)*1.25
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class.1, S.class),
                 mvrnorm(n=K/2, mu.class.2, S.class)))
  Sigmas <- abind(lapply(1:K, function(k) S.k), along=3)

  samp <- sim_gmm(mus=mus, Sigmas=Sigmas, n)
  return(list(X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=samp$Y,
              Z=c(rep(1, n/2), rep(2, n/2))))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc <- function(n, d, K=16, sigma=0) {
  # centers
  K.cent <- K/2
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)*1.25
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K.cent/2, mu.class.1, S.class),
                 mvrnorm(n=K.cent/2, mu.class.2, S.class)))

  ni <- n/K

  X <- do.call(rbind, lapply(1:K.cent, function(k) {
    X <- array(NaN, dim=c(ni*2, d))
    X[1:ni,] <- sweep(mgc.sims.2ball(ni, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni + 1):(2*ni),] <- sweep(mgc.sims.2sphere(ni, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X=X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=Y, Z=c(rep(1, n/2), rep(2, n/2))))
}
