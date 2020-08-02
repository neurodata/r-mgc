##----------------------------------------------
# Two Sample Tests
##----------------------------------------------
# should have current directory of this script as working directory
source('./shared_scripts.R')
no_cores = detectCores() - 1


## Two Sample Driver

test.two_sample <- function(X1, X2, Y, Z, dist.xfm=mgc.distance,
                            dist.params=list(method="euclidian"), dist.return=NULL,
                            remove.isolates=TRUE, nperm=100,
                            no_cores=1, alt="greater") {
  validated1 <- validator(X1, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)
  validated2 <- validator(X2, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)

  D1 <- validated1$D; X1 <- validated1$X; D2 <- validated2$D; X2 <- validated2$X; Y1 <- validated1$Y; Y2 <- validated2$Y; N <- nrow(D1)
  if (!all(Y1 == Y2)) {
    stop("The ids are not equal after removal of isolates from X1 and X2.")
  }
  if (nrow(D1) != nrow(D2) || nrow(D1) != length(Y1) || nrow(D1) != ncol(D1)) {
    stop("The distance matrices do not have the same number of elements.")
  }
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  } else if (no_cores >= 0.8*detectCores()) {
    warning("You have requested a number of cores near your machine's core count. Expected decreased performance.")
  }
  Xr1 <- lol.project.pca(X1, r=1)$Xr; Xr2 <- lol.project.pca(X2, r=1)$Xr
  # get observed difference in statistic of interest
  tr <- list(Discr=mgc:::discr.mnr(mgc:::discr.rdf(D1, Y1)) - mgc:::discr.mnr(mgc:::discr.rdf(D2, Y1)),
             PICC=icc.os(Xr1, Y1) - icc.os(Xr2, Y1),
             I2C2=i2c2.os(X1, Y1) - i2c2.os(X2, Y1),
             Kernel=ksamp.os(D1, Y1, is.dist=TRUE) - ksamp.os(D2, Y1, is.dist=TRUE),
             MMD=mmd.os(X1, Y1) - mmd.os(X2, Y1),
             HSIC=hsic.os(D1, Y1, is.dist=TRUE) - hsic.os(D2, Y1, is.dist=TRUE),
             DISCO=disco.os(D1, Y1, is.dist=TRUE) - disco.os(D2, Y1, is.dist=TRUE),
             FPI=fpi.os(X1, Y1, Z) - fpi.os(X2, Y1, Z))

  null.stats <- mclapply(1:nperm, function(i) {
    # generate null dataset for X1
    idx1 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda1 <- runif(N)
    Xn1 <- lambda1*X1[idx1[,1],] + (1 - lambda1)*X1[idx1[,2],]  # convex combination of elements of X1

    # generate null dataset for X2
    idx2 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda2 <- runif(N)
    Xn2 <- lambda2*X2[idx2[,1],] + (1 - lambda2)*X2[idx2[,2],]  # convex combination of elements of X2

    Xnr1 <- lol.project.pca(Xn1, r=1)$Xr; Xnr2 <- lol.project.pca(Xn2, r=1)$Xr
    DXn1 <- mgc.distance(Xn1); DXn2 <- mgc.distance(Xn2)
    # compute statistics of interest under the null
    D1.null <- discr.stat(DXn1, Y1, is.dist=TRUE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr
    D2.null <- discr.stat(DXn2, Y1, is.dist=TRUE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr
    ksamp1.null <- ksamp.os(DXn1, Y1, is.dist=TRUE)
    ksamp2.null <- ksamp.os(DXn2, Y1, is.dist=TRUE)
    disco1.null <- disco.os(DXn1, Y1, is.dist=TRUE)
    disco2.null <- disco.os(DXn2, Y1, is.dist=TRUE)
    hsic1.null <- disco.os(DXn1, Y1, is.dist=TRUE)
    hsic2.null <- disco.os(DXn2, Y1, is.dist=TRUE)
    mmd1.null <- mmd.os(Xn1, Y1)
    mmd2.null <- mmd.os(Xn2, Y2)
    icc1.null <- icc.os(Xnr1, Y1)
    icc2.null <- icc.os(Xnr2, Y1)
    i2c21.null <- i2c2.os(Xn1, Y1)
    i2c22.null <- i2c2.os(Xn2, Y1)
    fpi1.null <- fpi.os(Xn1, Y1, Z)
    fpi2.null <- fpi.os(Xn2, Y1, Z)
    return(list(Discr=c(D1.null - D2.null, D2.null - D1.null),
                PICC=c(icc1.null - icc2.null, icc2.null - icc1.null),
                I2C2=c(i2c21.null - i2c22.null, i2c22.null - i2c21.null),
                Kernel=c(ksamp1.null - ksamp2.null, ksamp2.null - ksamp1.null),
                MMD=c(mmd1.null - mmd2.null, mmd2.null - mmd1.null),
                DISCO=c(disco1.null - disco2.null, disco2.null - disco1.null),
                HSIC=c(hsic1.null - hsic2.null, hsic2.null - hsic1.null),
                FPI=c(fpi1.null - fpi2.null, fpi2.null - fpi1.null)))
  }, mc.cores=no_cores)

  # compute null distribution of difference between discriminabilities
  null.diff <- list(Discr=sapply(null.stats, function(x) x$Discr),
                    PICC=sapply(null.stats, function(x) x$PICC),
                    I2C2=sapply(null.stats, function(x) x$I2C2),
                    Kernel=sapply(null.stats, function(x) x$Kernel),
                    MMD=sapply(null.stats, function(x) x$MMD),
                    DISCO=sapply(null.stats, function(x) x$DISCO),
                    HSIC=sapply(null.stats, function(x) x$HSIC),
                    FPI=sapply(null.stats, function(x) x$FPI))

  return(do.call(rbind, lapply(names(tr), function(stat.name) {
    pval = mean(null.diff[[stat.name]] > tr[[stat.name]])*(nperm - 1)/nperm + 1/nperm
    if (is.na(pval)) {
      pval <- NaN
    }
    data.frame(stat.name=stat.name, stat=tr[[stat.name]],
               p.value=pval)
    }))
  )
}

## No Signal
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n=128, d=2, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp1 <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  samp2 <- sim_gmm_match(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), samp1$Y)
  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y,
              Z=sapply(unique(samp1$Y), function(y) return(1:sum(samp1$Y == y)))))
}

sim.parallel_rot_cigars <- function(n=128, d=2, sigma=0) {
  S.class <- diag(d)
  Sigma <- array(2, dim=c(d, d))
  diag(Sigma) <- 3

  mus <- cbind(c(0, 3), c(0, 0))
  samp1 <- sim_gmm(mus, Sigmas=abind(Sigma, Sigma, along=3), n)

  samp2 <- sim_gmm_match(mus, Sigmas=abind(Sigma, Sigma, along=3), samp1$Y)

  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y,
              Z=do.call(c, lapply(unique(samp1$Y), function(y) return(1:sum(samp1$Y == y))))))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable, and pipeline 1 is more discriminable
# than pipeline 2
# 2 classes
sim.linear_sig <- function(n, d, sigma=0) {
  S.class <- diag(d)
  Sigma <- diag(d)
  Sigma[1, 1] <- 2
  Sigma[-c(1), -c(1)] <- 1
  Sigma[1,1] <- 2
  mus.class <- mvrnorm(n=2, c(0,0), S.class)
  pi.k <- 0.5  # equal chance of a new sample being from class 1 or class 2
  samp1 <- sim_gmm(mus.class, Sigmas=abind(Sigma, Sigma, along=3), n)

  samp2 <- sim_gmm_match(mus.class, Sigmas=abind(Sigma, Sigma, along=3), samp1$Y)
  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y,
              Z=do.call(c, lapply(unique(samp1$Y), function(y) return(1:sum(samp1$Y == y))))))
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
  samp1 <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, n)
  samp2 <- sim_gmm_match(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, samp1$Y)

  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y,
              Z=do.call(c, lapply(unique(samp1$Y), function(y) return(1:sum(samp1$Y == y))))))
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
  X1 <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    # sample from crossed gaussians w p=0.5, 0.5 respectively
    sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, ni[k])
    # add individual-specific signal
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X2 <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rhos[k]
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rhos[k]
    sim <- sim_gmm(mus=cbind(rep(0, d), rep(0, d)), Sigmas=Sigmas, ni[k])
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  })) + array(rnorm(n*d), dim=c(n, d))*sigma

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))

  return(list(X1=X1, X2=X2, Y=Y,
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
  samp1 <- sim_gmm(mus=mus.class, Sigmas=Sigmas, n)
  samp2 <- sim_gmm_match(mus=mus.class, Sigmas=Sigmas, Ys=samp1$Y)
  return(list(X1=samp1$X, X2= samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y,
              Z=do.call(c, lapply(unique(samp1$Y), function(y) return(1:sum(samp1$Y == y))))))
}

sim.multiclass_ann_disc <- function(n, d, K=16, sigma=0) {
  # centers
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class, S.class)))

  ni <- rep(n/K, K)

  # individuals are either (1) a ball, or (2) a disc, around means
  X1 <- do.call(rbind, lapply(1:(K/2), function(k) {
    n.ball <- ni[2*(k-1)+1]; n.disc <- ni[2*k]
    X <- array(NaN, dim=c((n.ball + n.disc), d))
    X[1:n.ball,] <- sweep(mgc.sims.2ball(n.ball, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(n.ball + 1):(n.ball + n.disc),] <- sweep(mgc.sims.2sphere(n.disc, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))
  X2 <- do.call(rbind, lapply(1:(K/2), function(k) {
    n.ball <- ni[2*(k-1)+1]; n.disc <- ni[2*k]
    X <- array(NaN, dim=c((n.ball + n.disc), d))
    X[1:n.ball,] <- sweep(mgc.sims.2ball(n.ball, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(n.ball + 1):(n.ball + n.disc),] <- sweep(mgc.sims.2sphere(n.disc, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))
  return(list(X1=X1, X2=X2 + array(rnorm(n*d), dim=c(n, d))*sigma, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
}

# 8 pairs of annulus/discs
sim.multiclass_ann_disc2 <- function(n, d, sigma=0) {

  mus <- cbind(c(0, 0))

  ni <- rep(n/2, 2)

  X <- array(NaN, dim=c(n, d))
  X[1:ni[1],] <- sweep(mgc.sims.2ball(ni[1], d, r=1, cov.scale=0.1), 2, mus[,1], "+")
  X[(ni[1] + 1):n,] <- sweep(mgc.sims.2sphere(ni[2], r=1.5, d=d, cov.scale=0.1), 2, mus[,1], "+")

  Y <- c(rep(1, ni[1]), rep(2, ni[2]))

  X2 <- array(NaN, dim=c(n, d))
  X2[1:ni[1],] <- sweep(mgc.sims.2ball(ni[1], d, r=1, cov.scale=0.1), 2, mus[,1], "+")
  X2[(ni[1] + 1):n,] <- sweep(mgc.sims.2sphere(ni[2], r=1.5, d=d, cov.scale=0.1), 2, mus[,1], "+")

  return(list(X1=X, X2=X2 + array(rnorm(n*d), dim=c(n, d))*sigma, Y=Y,
              Z=do.call(c, lapply(unique(Y), function(y) return(1:sum(Y == y))))))
}

sim.xor2 <- function(n, d, sigma=0) {
  mus <- cbind(c(0, 0), c(1,1), c(1, 0), c(0, 1))
  Y <- rep(1:ncol(mus), n/ncol(mus))
  X1 <- mvrnorm(n=n, mu=c(0, 0), Sigma=.1*diag(d)) + t(mus[,Y])
  X2 <- mvrnorm(n=n, mu=c(0, 0), Sigma=.1*diag(d)) + t(mus[,Y])
  Y <- floor((Y-1)/2)
  Z <- Y
  Z[Y == 0] <- sample(1:(n/2), replace=FALSE, size=n/2)
  Z[Y == 1] <- sample(1:(n/2), replace=FALSE, size=n/2)
  return(list(X1=X1, X2=X2 + array(rnorm(n*d), dim=c(n, d))*sigma, Y=Y,
              Z=Z))
}

## --------------------------------------
# Driver
## --------------------------------------
n <- 128; d <- 2
nrep <- 200
n.sigma <- 15

simulations <- list(sim.no_signal, sim.crossed_sig2,
                    sim.multiclass_gaussian, sim.multiclass_ann_disc2, sim.xor2)
sims.sig.max <- c(10, 1, 1, 1, .2)
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

list.results.ts <- mclapply(1:length(experiments), function(i) {
  exper <- experiments[[i]]
  print(i)
  sim <- simpleError("Fake Error"); att = 0
  while(inherits(sim, "error") && att <= 50) {
    sim <- tryCatch({
      simu <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
      if (dim(simu$X1)[1] != length(simu$Y) || dim(simu$X2)[1] != length(simu$Y)) {
        stop("Error")
      }
      simu
    }, error=function(e) e)
    att <- att + 1
  }
  res <- test.two_sample(sim$X1, sim$X2, sim$Y, sim$Z, nperm=100)
  res$sim.name <- exper$sim.name; res$n <- n; res$d <- d; res$i <- exper$i
  res$sigma <- exper$sigma
  return(res)
}, mc.cores=no_cores)

ts.results <- do.call(rbind, list.results.ts)
saveRDS(list(ts.results=ts.results, list.results=list.results.ts), '../data/sims/discr_sims_ts.rds')

