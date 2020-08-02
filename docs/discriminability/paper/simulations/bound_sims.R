##----------------------------------------------
# Bound Simulation
##----------------------------------------------
# should have current directory of this script as working directory
source('./shared_scripts.R')
no_cores = detectCores() - 1

compute_bayes <- function(X, Z, npts=100) {
  # compute the range of the dimensions occupied across
  # both classes
  range.x <- c(min(X[,1]), max(X[,1]))
  range.y <- c(min(X[,2]), max(X[,2]))

  # compute 2d density estimate
  gr1.dens=kde2d(X[Z == 1,1], X[Z == 1,2], lims=c(range.x, range.y), n=npts)
  gr1.dens$z <- gr1.dens$z/sum(gr1.dens$z)  # normalize so it is a pdf
  gr2.dens=kde2d(X[Z == 2,1], X[Z == 2,2], lims=c(range.x, range.y), n=npts)
  gr2.dens$z <- gr2.dens$z/sum(gr2.dens$z)

  # compute the density we would get wrong by using the most likely
  # points
  bayes <- sum(pmin(gr1.dens$z, gr2.dens$z))*0.5
  return(bayes)
}

## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n=128, d=2, n.bayes=10000, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  samp$X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma

  samp.bayes <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n.bayes)
  samp.bayes$X=samp.bayes$X + array(rnorm(n.bayes*d), dim=c(n.bayes, d))*sigma

  DX <- mgc.distance(samp$X)
  Z <- do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))
  return(list(Discr=discr.stat(DX, samp$Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(samp$X, r=1)$Xr, samp$Y),
              I2C2=i2c2.os(samp$X, samp$Y), Kernel=ksamp.os(DX, samp$Y, is.dist=TRUE), FPI=fpi.os(samp$X, samp$Y, Z),
              MMD=mmd.os(samp$X, samp$Y), HSIC=hsic.os(DX, samp$Y, is.dist=TRUE), DISCO=disco.os(DX, samp$Y, is.dist=TRUE),
              bayes=compute_bayes(samp.bayes$X, samp.bayes$Y)))
}

sim.parallel_rot_cigars <- function(n=128, d=2, n.bayes=10000, n.pts=100, sigma=0) {
  S.class <- diag(d)
  Sigma <- array(2, dim=c(d, d))
  diag(Sigma) <- 3

  mus <- cbind(c(0, 3), c(0, 0))
  samp <- sim_gmm(mus, Sigmas=abind(Sigma, Sigma, along=3), n)
  samp$X <- samp$X + array(rnorm(n*d), dim=c(n, d))*sigma

  samp.bayes <- sim_gmm(mus, Sigmas=abind(Sigma, Sigma, along=3), n.bayes)
  samp.bayes$X <- samp.bayes$X + array(rnorm(n.bayes*d), dim=c(n.bayes, d))*sigma

  DX <- mgc.distance(samp$X)
  Z <- do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))
  return(list(Discr=discr.stat(DX, samp$Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(samp$X, r=1)$Xr, samp$Y),
              I2C2=i2c2.os(samp$X, samp$Y), Kernel=ksamp.os(DX, samp$Y, is.dist=TRUE), FPI=fpi.os(samp$X, samp$Y, Z),
              MMD=mmd.os(samp$X, samp$Y), HSIC=hsic.os(DX, samp$Y, is.dist=TRUE), DISCO=disco.os(DX, samp$Y, is.dist=TRUE),
              bayes=compute_bayes(samp.bayes$X, samp.bayes$Y)))
}


## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n=128, d=2, n.bayes=10000, n.pts=100, sigma=0) {
  S.class <- diag(d)
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus.class <- mvrnorm(n=2, c(0,0), S.class)
  samp <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n)
  samp$X <- samp$X + array(rnorm(n*d), dim=c(n, d))*sigma

  samp.bayes <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n.bayes)
  samp.bayes$X <- samp.bayes$X + array(rnorm(n.bayes*d), dim=c(n.bayes, d))*sigma

  DX <- mgc.distance(samp$X)
  Z <- do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))
  return(list(Discr=discr.stat(DX, samp$Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(samp$X, r=1)$Xr, samp$Y),
              I2C2=i2c2.os(samp$X, samp$Y), Kernel=ksamp.os(DX, samp$Y, is.dist=TRUE), FPI=fpi.os(samp$X, samp$Y, Z),
              MMD=mmd.os(samp$X, samp$Y), HSIC=hsic.os(DX, samp$Y, is.dist=TRUE), DISCO=disco.os(DX, samp$Y, is.dist=TRUE),
              bayes=compute_bayes(samp.bayes$X, samp.bayes$Y)))
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

  ni <- n/K
  rhos <- runif(K, min=-.2, max=.2)
  X <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, ni)
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))

  # probability of being each individual is 1/K
  ni.bayes <- n.bayes/K

  X.bayes <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, ni.bayes[k])
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X.bayes <- X.bayes + array(rnorm(n.bayes*d)*sigma, dim=c(n.bayes, d))

  Y.bayes <- do.call(c, lapply(1:K, function(k) rep(k, ni.bayes)))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- as.numeric(sapply(1:K, function(k) {
    return(mus.class[1,k] <= 0)
  })) + 1
  Z.bayes <- do.call(c, lapply(1:(K/2), function(k) {
    return(rep(mus.z[k], 2*ni.bayes))
  }))
  DX <- mgc.distance(X)
  Z <- do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))
  return(list(Discr=discr.stat(DX, Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(X, r=1)$Xr, Y),
              I2C2=i2c2.os(X, Y), Kernel=ksamp.os(DX, Y, is.dist=TRUE), FPI=fpi.os(X, Y, Z),
              MMD=mmd.os(X, Y), HSIC=hsic.os(DX, Y, is.dist=TRUE), DISCO=disco.os(DX, Y, is.dist=TRUE),
              bayes=compute_bayes(X.bayes, Z.bayes)))
}

## Crossed Signal Difference
# a simulation where classes are crossed but distinguishable
# also contains correlation btwn dimensions
# 2 classes
sim.crossed_sig2 <- function(n=128, d=2, n.bayes=10000, sigma=0) {
  # class mus
  K=2
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus.class <- t(mvrnorm(n=K, mu.class, S.class))

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))  # covariances are orthogonal
  mus=cbind(rep(0, d), rep(0, d))

  rho <- runif(1, min=-.2, max=.2)

  # add random correlation
  Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
  Sigmas[1,2,1] <- Sigmas[2,1,1] <- rho
  Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rho
  # sample from crossed gaussians w p=0.5, 0.5 respectively
  sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, n)

  X <- sim$X + array(rnorm(n*d)*sigma, dim=c(n, d))
  Y <- sim$Y

  # Bayes Sim
  sim.bayes <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, n.bayes)
  X.bayes <- sim.bayes$X + array(rnorm(n.bayes*d)*sigma, dim=c(n.bayes, d))
  Y.bayes <- sim.bayes$Y

  DX <- mgc.distance(X)
  Z <- do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))
  return(list(Discr=discr.stat(DX, Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(X, r=1)$Xr, Y),
              I2C2=i2c2.os(X, Y), Kernel=ksamp.os(DX, Y, is.dist=TRUE), FPI=fpi.os(X, Y, Z),
              MMD=mmd.os(X, Y), HSIC=hsic.os(DX, Y, is.dist=TRUE), DISCO=disco.os(DX, Y, is.dist=TRUE),
              bayes=compute_bayes(X.bayes, Y.bayes)))
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
  samp <- sim_gmm(mus=mus.class, Sigmas=Sigmas, n)
  samp$X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d))

  # sample individuals w.p. 1/K
  samp.bayes <- sim_gmm(mus=mus.class, Sigmas=Sigmas, n.bayes)
  samp.bayes$X=samp.bayes$X + array(rnorm(n.bayes*d)*sigma, dim=c(n.bayes, d))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- as.numeric(sapply(1:K, function(k) {
    return(mus.class[1,k] <= 0)
  })) + 1
  Z.bayes <- mus.z[samp.bayes$Y]

  DX <- mgc.distance(samp$X)
  Z <- do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))
  return(list(Discr=discr.stat(DX, samp$Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(samp$X, r=1)$Xr, samp$Y),
              I2C2=i2c2.os(samp$X, samp$Y), Kernel=ksamp.os(DX, samp$Y, is.dist=TRUE), FPI=fpi.os(samp$X, samp$Y, Z),
              MMD=NaN, HSIC=hsic.os(DX, samp$Y, is.dist=TRUE), DISCO=disco.os(DX, samp$Y, is.dist=TRUE),
              bayes=compute_bayes(samp.bayes$X, Z.bayes)))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc <- function(n, d, K=16, n.bayes=10000, sigma=0) {
  # centers
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class, S.class)))

  # probability of being each individual is 1/K
  ni <- n/K

  # individuals are either (1) a ball, or (2) a disc, around means
  X <- do.call(rbind, lapply(1:(K/2), function(k) {
    X <- array(NaN, dim=c((2*ni), d))
    X[1:ni,] <- sweep(mgc.sims.2ball(ni, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni + 1):(2*ni),] <- sweep(mgc.sims.2sphere(ni, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))
  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))

  # probability of being each individual is 1/K
  ni.bayes <- n.bayes/K

  # individuals are either (1) a ball, or (2) a disc, around means
  X.bayes <- do.call(rbind, lapply(1:(K/2), function(k) {
    X <- array(NaN, dim=c((2*ni.bayes), d))
    X[1:ni.bayes,] <- sweep(mgc.sims.2ball(ni.bayes, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni.bayes + 1):(2*ni.bayes),] <- sweep(mgc.sims.2sphere(ni.bayes, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))
  Y.bayes <- do.call(c, lapply(1:K, function(k) rep(k, ni.bayes)))
  # assign cluster centers based on <= 0 in first dimension
  mus.z <- as.numeric(sapply(1:(K/2), function(k) {
    return(mus[1,k] <= 0)
  })) + 1
  Z.bayes <- do.call(c, lapply(1:(K/2), function(k) {
    return(rep(mus.z[k], 2*ni.bayes))
  }))

  DX <- mgc.distance(X)
  Z <- do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))
  return(list(Discr=discr.stat(DX, Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(X, r=1)$Xr, Y),
              I2C2=i2c2.os(X, Y), Kernel=ksamp.os(DX, Y, is.dist=TRUE), FPI=fpi.os(X, Y, Z),
              MMD=NaN, HSIC=hsic.os(DX, Y, is.dist=TRUE), DISCO=disco.os(DX, Y, is.dist=TRUE),
              bayes=compute_bayes(X.bayes, Z.bayes)))
}


# 8 pairs of annulus/discs
sim.multiclass_ann_disc2 <- function(n, d, n.bayes=5000, sigma=0) {

  mus <- cbind(c(0, 0))

  # probability of being each individual is 1/K
  ni <- rep(n/2, 2)

  X <- array(NaN, dim=c(n, d))
  X[1:ni[1],] <- sweep(mgc.sims.2ball(ni[1], d, r=1, cov.scale=0.1), 2, mus[,1], "+")
  X[(ni[1] + 1):n,] <- sweep(mgc.sims.2sphere(ni[2], r=1.5, d=d, cov.scale=0.1), 2, mus[,1], "+")

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- c(rep(1, ni[1]), rep(2, ni[2]))

  # probability of being each individual is 1/K
  ni.bayes <- rep(n.bayes/2, 2)

  X.bayes <- array(NaN, dim=c(n.bayes, d))
  X.bayes[1:ni.bayes[1],] <- sweep(mgc.sims.2ball(ni.bayes[1], d, r=1, cov.scale=0.1),
                                   2, mus[,1], "+")
  X.bayes[(ni.bayes[1] + 1):n.bayes,] <- sweep(mgc.sims.2sphere(ni.bayes[2], r=1, d=d, cov.scale=0.1),
                                               2, mus[,1], "+")
  X.bayes <- X.bayes + array(rnorm(n*d)*sigma, dim=c(n.bayes, d))
  Y.bayes <- c(rep(1, ni.bayes[1]), rep(2, ni.bayes[2]))

  DX <- mgc.distance(X)
  Z <- do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))
  return(list(Discr=discr.stat(DX, Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(X, r=1)$Xr, Y),
              I2C2=i2c2.os(X, Y), Kernel=ksamp.os(DX, Y, is.dist=TRUE), FPI=fpi.os(X, Y, Z),
              MMD=NaN, HSIC=hsic.os(DX, Y, is.dist=TRUE), DISCO=disco.os(DX, Y, is.dist=TRUE),
              bayes=compute_bayes(X.bayes, Y.bayes)))
}


sim.xor2 <- function(n, d, n.bayes=10000, sigma=0) {
  mus <- cbind(c(0, 0), c(1,1), c(1, 0), c(0, 1))
  Y <- rep(1:ncol(mus), n/ncol(mus))
  X <- mvrnorm(n=n, mu=c(0, 0), Sigma=sigma*diag(d)) + t(mus[,Y])
  Y <- floor((Y-1)/2) + 1

  Y.bayes <- rep(1:ncol(mus), n.bayes/ncol(mus))
  X.bayes <- mvrnorm(n=n.bayes, mu=c(0, 0), Sigma=sigma*diag(d)) + t(mus[,Y.bayes])
  Y.bayes <- floor((Y.bayes-1)/2) + 1

  DX <- mgc.distance(X)
  Z <- Y
  Z[Y == 0] <- sample(1:(n/2), replace=FALSE, size=n/2)
  Z[Y == 1] <- sample(1:(n/2), replace=FALSE, size=n/2)
  return(list(Discr=discr.stat(DX, Y, is.dist=TRUE)$discr, PICC=icc.os(lol.project.pca(X, r=1)$Xr, Y),
              I2C2=i2c2.os(X, Y), Kernel=ksamp.os(DX, Y, is.dist=TRUE), FPI=fpi.os(X, Y, Z),
              MMD=mmd.os(X, Y), HSIC=hsic.os(DX, Y, is.dist=TRUE), DISCO=disco.os(DX, Y, is.dist=TRUE),
              bayes=compute_bayes(X.bayes, Y.bayes)))
}

n <- 128; d <- 2
nrep <- 200
# 200, 15
n.sigma <- 15

simulations <- list(sim.no_signal, sim.crossed_sig2,
                    sim.multiclass_gaussian, sim.multiclass_ann_disc2, sim.xor2)
sims.sig.max <- c(10, 2, 2, 1, 0.5)
sims.sig.min <- c(0, 0, 0, 0, 0)
names(simulations) <- names(sims.sig.max) <- names(sims.sig.min) <-
  c("No Signal", "Cross", "Gaussian", "Ball/Circle", "XOR")

experiments <- do.call(c, lapply(names(simulations), function(sim.name) {
  do.call(c, lapply(seq(from=sims.sig.min[sim.name], to=sims.sig.max[sim.name],
                        length.out=n.sigma), function(sigma) {
                          lapply(1:nrep, function(i) {
                            return(list(i=i, sim=simulations[[sim.name]], sim.name=sim.name, sigma=sigma))
                          })
                        }))
}))

list.results.bound <- mclapply(1:length(experiments), function(i) {
  exper <- experiments[[i]]
  print(i)
  sim <- simpleError("Fake Error"); att = 0
  while(inherits(sim, "error") && att <= 50) {
    sim <- tryCatch({
      simu <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
      simu
    }, error=function(e) {e})
    att <- att + 1
  }
  res <- data.frame(sim.name=exper$sim.name, n=n, d=d, i=exper$i, algorithm=names(sim),
                    sigma=exper$sigma, value=as.numeric(do.call(c, sim)))
  return(res)
}, mc.cores=no_cores)

bound.results <- do.call(rbind, list.results.bound)
saveRDS(list(bound.results=bound.results, list.results=list.results.bound), '../data/sims/discr_sims_bound.rds')

