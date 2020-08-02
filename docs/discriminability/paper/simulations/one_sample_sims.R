##----------------------------------------------
# One Sample Tests
##----------------------------------------------
# should have current directory of this script as working directory
source('./shared_scripts.R')
no_cores = detectCores() - 1

## One Sample Driver
test.one_sample <- function(X, Y, Z, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                            dist.return=NULL, remove.isolates=TRUE, nperm=100, no_cores=1) {

  validated <- validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                     remove.isolates=remove.isolates)

  D <- validated$D; Y <- validated$Y; X <- validated$X; N <- nrow(D)
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  } else if (no_cores >= 0.8*detectCores()) {
    warning("You have requested a number of cores near your machine's core count. Expected decreased performance.")
  }
  Xr <- lol.project.pca(X, 1)$Xr
  # compute references for the statistics
  tr <- list(Discr=discr.stat(D, Y, is.dist=TRUE)$discr,
             PICC=icc.os(Xr, Y),
             I2C2=i2c2.os(X, Y),
             Kernel=ksamp.os(D, Y, is.dist=TRUE),
             MMD=mmd.os(X, Y),
             HSIC=hsic.os(D, Y, is.dist=TRUE),
             DISCO=disco.os(D, Y, is.dist=TRUE),
             FPI=fpi.os(X, Y, Z))

  nr <- mclapply(1:nperm, function(i) {
    perm.idx <- sample(N)
    perm.Y <- Y[perm.idx]
    perm.Z <- Z[perm.idx]
    return(list(Discr=discr.stat(D, perm.Y, is.dist=TRUE)$discr,
                PICC=icc.os(Xr, perm.Y),
                I2C2=i2c2.os(X, perm.Y),
                Kernel=ksamp.os(D, perm.Y, is.dist=TRUE),
                MMD=mmd.os(X, perm.Y),
                HSIC=hsic.os(D, perm.Y, is.dist=TRUE),
                DISCO=disco.os(D, perm.Y,is.dist=TRUE),
                FPI=fpi.os(X, perm.Y, perm.Z)))
  }, mc.cores=no_cores)

  null.stats <- list(Discr=lapply(nr, function(x) x$Discr),
                     PICC=lapply(nr, function(x) x$PICC),
                     I2C2=lapply(nr, function(x) x$I2C2),
                     Kernel=lapply(nr, function(x) x$Kernel),
                     MMD=lapply(nr, function(x) x$MMD),
                     HSIC=lapply(nr, function(x) x$HSIC),
                     DISCO=lapply(nr, function(x) x$DISCO),
                     FPI=lapply(nr, function(x) x$FPI))

  return(do.call(rbind, lapply(names(tr), function(stat.name) {
    pval = mean(null.stats[[stat.name]] > tr[[stat.name]])*(nperm - 1)/nperm + 1/nperm
    if (is.na(pval)) {
      pval <- NaN
    }
    data.frame(stat.name=stat.name, stat=tr[[stat.name]],
               p.value=pval)
  })))
}

## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n=128, d=2, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  samp$X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma
  return(list(X=samp$X, Y=samp$Y,
              Z=do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))))
}

sim.parallel_rot_cigars <- function(n=128, d=2, sigma=0) {
  S.class <- diag(d)
  Sigma <- array(2, dim=c(d, d))
  diag(Sigma) <- 3

  mus <- cbind(c(0, 2.5), c(0, 0))
  samp <- sim_gmm(mus, Sigmas=abind(Sigma, Sigma, along=3), n)

  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y,
              Z=do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n, d, sigma=0) {
  S.class <- diag(d)
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus.class <- mvrnorm(n=2, c(0,0), S.class)
  pi.k <- 0.5  # equal chance of a new sample being from class 1 or class 2
  samp <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n)
  samp$X <- samp$X + array(rnorm(n*d), dim=c(n, d))*sigma
  return(list(X=samp$X, Y=samp$Y,
              Z=do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))))
}

sim.crossed_sig2 <- function(n=128, d=2, sigma=0) {
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
  return(list(X=X, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
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

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))  # covariances are orthogonal
  mus=cbind(rep(0, d), rep(0, d))

  ni <- rep(n/K, K)
  rhos <- runif(K, min=-.2, max=.2)
  X <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, ni[k])
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
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
  samp <- sim_gmm(mus=mus.class, Sigmas=Sigmas, n)
  samp$X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d))
  return(list(X=samp$X, Y=samp$Y,
              Z=do.call(c, lapply(unique(samp$Y), function(y) return(1:sum(samp$Y == y))))))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc <- function(n, d, K=16, sigma=0) {
  # centers
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class, S.class)))

  ni <- rep(n/K, K)

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
  return(list(X=X, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc2 <- function(n, d, sigma=0) {

  mus <- cbind(c(0, 0))

  ni <- rep(n/2, 2)

  X <- array(NaN, dim=c(n, d))
  X[1:ni[1],] <- sweep(mgc.sims.2ball(ni[1], d, r=1, cov.scale=0.1), 2, mus[,1], "+")
  X[(ni[1] + 1):n,] <- sweep(mgc.sims.2sphere(ni[2], r=1.5, d=d, cov.scale=0.1), 2, mus[,1], "+")

  X <- X + array(rnorm(n*d)*sigma, dim=c(n, d))

  Y <- c(rep(1, ni[1]), rep(2, ni[2]))
  return(list(X=X, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
}

sim.xor2 <- function(n, d, sigma=0) {

  mus <- cbind(c(0, 0), c(1,1), c(1, 0), c(0, 1))
  Y <- rep(1:ncol(mus), n/ncol(mus))
  X <- mvrnorm(n=n, mu=c(0, 0), Sigma=sigma*diag(d)) + t(mus[,Y])
  Y <- floor((Y-1)/2)
  Z <- Y
  Z[Y == 0] <- sample(1:(n/2), replace=FALSE, size=n/2)
  Z[Y == 1] <- sample(1:(n/2), replace=FALSE, size=n/2)
  return(list(X=X, Y=Y, Z=Z))
}
## -------------------------
# Driver
## -------------------------
n <- 128; d <- 2
nrep <- 200 # 500
n.sigma <- 15 # 15
nperm <- 100

simulations <- list(sim.no_signal, sim.crossed_sig2,
                    sim.multiclass_gaussian, sim.multiclass_ann_disc2,
                    sim.xor2)
sims.sig.max <- c(20, 2, 20, 2, 0.8)
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


list.results.os <- mclapply(1:length(experiments), function(i) {
  exper <- experiments[[i]]
  print(i)
  sim <- simpleError("Fake Error"); att = 0
  while(inherits(sim, "error") && att <= 50) {
    sim <- tryCatch({
      simu <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
      if (dim(simu$X)[1] != length(simu$Y)) {
        stop("Error")
      }
      simu
    }, error=function(e) print(e))
    att <- att + 1
  }
  res <- test.one_sample(sim$X, sim$Y, sim$Z, nperm=nperm)
  res$sim.name <- exper$sim.name; res$n <- n; res$d <- d; res$i <- exper$i
  res$sigma <- exper$sigma
  return(res)
}, mc.cores=no_cores)

bound.results.os <- do.call(rbind, list.results.os)
saveRDS(list(os.results=bound.results.os, list.results=list.results.os), '../data/sims/discr_sims_os.rds')

