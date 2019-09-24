##----------------------------------------------
# Bound Simulation
##----------------------------------------------

require(MASS)
require(mvtnorm)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(lolR)
require(abind)
require(emdbook)
no_cores = detectCores() - 1

## -------------------------------------------
# Algorithms
## -------------------------------------------
# one-way ICC
icc.os <- function(x, y) {
  data <- data.frame(x=x, y=y)
  fit <- anova(aov(x ~ y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  a <- length(unique(y))
  tmp.outj <- as.numeric(aggregate(x ~ y, data=data, FUN = length)$x)
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  return(r)
}

# I2C2 wrapper
i2c2.os <- function(X, Y) {
  return(I2C2.original(y=X, id=Y, visit=rep(1, length(Y)), twoway=FALSE)$lambda)
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


compute_bayes <- function(X, Z, npts=100) {
  # compute the range of the dimensions occupied across
  # both classes
  range.x <- c(min(samp$X[,1]), max(samp$X[,1]))
  range.y <- c(min(samp$X[,2]), max(samp$X[,2]))

  # compute 2d density estimate
  gr1.dens=kde2d(samp$X[samp$Z == 1,1], samp$X[samp$Z == 1,2], lims=c(range.x, range.y), n=npts)
  gr1.dens$z <- test.first$z/sum(test.first$z)  # normalize so it is a pdf
  gr2.dens=kde2d(samp$X[samp$Z == 2,1], samp$X[samp$Z == 2,2], lims=c(range.x, range.y), n=npts)
  test.sec$z <- test.sec$z/sum(test.sec$z)

  gr2.dens <- array(NaN, dim=c(npts, npts))
  # compute the density we would get wrong by using the most likely
  # points
  bayes <- sum(pmin(test.first$z, test.sec$z))*0.5
  return(bayes)
}

## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n, priors=NULL) {
  K <- dim(mus)[2]
  if (is.null(priors)) {
    priors <- rep(1/K, K)
  }
  ni <- rowSums(rmultinom(n, 1, prob=priors))
  X <- do.call(rbind, lapply(1:K, function(k) mvrnorm(n=ni[k], mus[,k], Sigmas[,,k])))
  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X, Y=Y))
}

## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n=128, d=2, n.bayes=10000, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n, priors=c(0.5,0.5))
  samp$X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma

  samp.bayes <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n.bayes, priors=c(0.5,0.5))
  samp.bayes$X=samp.bayes$X + array(rnorm(n.bayes*d), dim=c(n.bayes, d))*sigma

  return(list(discr=discr.stat(samp$X, samp$Y)$discr, icc=icc.os(samp$X, samp$Y),
              i2c2=i2c2.os(samp$X, samp$Y), bayes=compute_bayes(samp.bayes$X, samp.bayes$Y)))
}


## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n=128, d=2, n.bayes=10000, n.pts=100,sigma=0) {
  S.class <- diag(d)
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus.class <- mvrnorm(n=2, c(0,0), S.class)
  pi.k <- 0.5  # equal chance of a new sample being from class 1 or class 2
  samp <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n, priors=c(pi.k, pi.k))
  samp$X <- samp$X + array(rnorm(n*d), dim=c(n, d))*sigma

  samp.bayes <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n.bayes, priors=c(pi.k, pi.k))
  samp.bayes$X <- samp.bayes$X + array(rnorm(n.bayes*d), dim=c(n.bayes, d))*sigma

  return(list(discr=discr.stat(samp$X, samp$Y)$discr, icc=icc.os(samp$X, samp$Y),
              i2c2=i2c2.os(samp$X, samp$Y), bayes=compute_bayes(samp.bayes$X, samp.bayes$Y)))
}

## Crossed Signal Difference
# a simulation where classes are crossed but distinguishable
# also contains correlation btwn dimensions
# 2 classes
sim.crossed_sig <- function(n=128, d=2, K=16, n.bayes=10000, sigma=0) {
  # class mus
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus.class <- t(mvrnorm(n=K, mu.class, S.class))

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))  # covariances are orthogonal
  mus=cbind(rep(0, d), rep(0, d))

  # probability of being each individual is 1/K
  ni <- rowSums(rmultinom(n, 1, prob=rep(1/K, K)))
  rhos <- runif(K, min=-1, max=1)
  X <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- mgc:::mgc.sims.sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas,
                                  ni[k], priors=c(0.5, 0.5))
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))

  # probability of being each individual is 1/K
  ni.bayes <- rowSums(rmultinom(n.bayes, 1, prob=rep(1/K, K)))

  X.bayes <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    *sqrt(sum(diag(Sigma.1)))
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- mgc:::mgc.sims.sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas,
                                  ni.bayes[k], priors=c(0.5, 0.5))
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X.bayes <- X.bayes + array(rnorm(n.bayes*d)*sigma, dim=c(n.bayes, d))

  Y.bayes <- do.call(c, lapply(1:K, function(k) rep(k, ni.bayes[k])))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- sapply(1:(K/2), function(k) {
    return(mus[1,k] <= 0)
  })
  Z.bayes <- do.call(c, lapply(1:(K/2), function(k) {
    return(rep(mus.z[k], ni.bayes[2*k - 1] + ni.bayes[2*k]))
  }))
  return(list(discr=discr.stat(samp$X, samp$Y)$discr, icc=icc.os(samp$X, samp$Y),
              i2c2=i2c2.os(samp$X, samp$Y), bayes=compute_bayes(samp.bayes$X, Z.bayes)))
}

## Samples from Multiclass Gaussians
# a simulation where there are multiple classes present, and a correlation structure
# 2 classes
sim.multiclass_gaussian <- function(n, d, K=16, n.bayes=10000, sigma=0) {
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
  samp$X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d))

  # sample individuals w.p. 1/K
  samp.bayes <- mgc:::mgc.sims.sim_gmm(mus=mus, Sigmas=Sigmas, n.bayes, priors=rep(1/K, K))
  samp.bayes$X=samp.bayes$X + array(rnorm(n.bayes*d)*sigma, dim=c(n.bayes, d))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- sapply(1:(K/2), function(k) {
    return(mus[1,k] <= 0)
  })
  Z.bayes <- do.call(c, lapply(1:(K/2), function(k) {
    return(rep(mus.z[k], ni.bayes[2*k - 1] + ni.bayes[2*k]))
  }))

  return(list(discr=discr.stat(samp$X, samp$Y)$discr, icc=icc.os(samp$X, samp$Y),
              i2c2=i2c2.os(samp$X, samp$Y), bayes=compute_bayes(samp.bayes$X, Z.bayes)))
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
  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))

  # probability of being each individual is 1/K
  ni.bayes <- rowSums(rmultinom(n.bayes, 1, prob=rep(1/K, K)))

  # individuals are either (1) a ball, or (2) a disc, around means
  X <- do.call(rbind, lapply(1:(K/2), function(k) {
    n.ball <- ni.bayes[2*(k-1)+1]; n.disc <- ni.bayes[2*k]
    X <- array(NaN, dim=c((n.ball + n.disc), d))
    X[1:n.ball,] <- sweep(mgc.sims.2ball(n.ball, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(n.ball + 1):(n.ball + n.disc),] <- sweep(mgc.sims.2sphere(n.disc, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))
  Y.bayes <- do.call(c, lapply(1:K, function(k) rep(k, ni.bayes[k])))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- sapply(1:(K/2), function(k) {
    return(mus[1,k] <= 0)
  })
  Z.bayes <- do.call(c, lapply(1:(K/2), function(k) {
    return(rep(mus.z[k], ni.bayes[2*k - 1] + ni.bayes[2*k]))
  }))

  return(list(discr=discr.stat(samp$X, samp$Y)$discr, icc=icc.os(samp$X, samp$Y),
              i2c2=i2c2.os(samp$X, samp$Y), bayes=compute_bayes(samp.bayes$X, Z.bayes)))
}
n <- 128; d <- 2
nrep <- 200
n.sigma <- 15

simulations <- list(sim.no_signal, sim.linear_sig, sim.crossed_sig,
                    sim.multiclass_gaussian, sim.multiclass_ann_disc)
sims.sig.max <- c(10, 2, 2, 2, 1)
sims.sig.min <- c(0, 0, 0, 0, 0)
names(simulations) <- names(sims.sig.max) <- names(sims.sig.min) <-
  c("No Signal", "Linear", "Cross", "Gaussian", "Annulus/Disc")

experiments <- do.call(c, lapply(names(simulations), function(sim.name) {
  do.call(c, lapply(seq(from=sims.sig.min[sim.name], to=sims.sig.max[sim.name],
                        length.out=n.sigma), function(sigma) {
                          lapply(1:nrep, function(i) {
                            return(list(i=i, sim=simulations[[sim.name]], sim.name=sim.name, sigma=sigma))
                          })
                        }))
}))

list.results.bound <- mclapply(experiments, function(exper) {
  er = 0
  while(er <= 50) {
    tryCatch({
      sim <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
      er <- 51
    }, error=function(e) {er <- er + 1})
  }
  res <- data.frame(sim.name=exper$sim.name, n=n, d=d, i=exper$i, algorithm=names(sim),
                    sigma=exper$sigma, value=as.numeric(do.call(c, sim)))
  return(res)
}, mc.cores=no_cores)

bound.results <- do.call(rbind, list.results.ts)
saveRDS(list(ts.results=ts.results, list.results=list.results.ts), '../data/sims/discr_sims_ts.rds')

