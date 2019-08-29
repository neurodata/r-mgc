##----------------------------------------------
# Two Sample Tests
##----------------------------------------------

require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(lolR)
require(abind)
no_cores = detectCores() - 1

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

## Two Sample Driver

test.two_sample <- function(X1, X2, Y, dist.xfm=mgc.distance,
                            dist.params=list(method="euclidian"), dist.return=NULL,
                            remove.isolates=TRUE, nperm=100,
                            no_cores=1, alt="greater") {
  validated1 <- mgc:::discr.validator(X1, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)
  validated2 <- mgc:::discr.validator(X2, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)

  D1 <- validated1$D; D2 <- validated2$D; Y1 <- validated1$Y; Y2 <- validated2$Y; N <- nrow(D1)
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
  tr <- list(discr=mgc:::discr.mnr(mgc:::discr.rdf(D1, Y1)) - mgc:::discr.mnr(mgc:::discr.rdf(D2, Y1)),
             icc=icc.os(Xr1, Y1) - icc.os(Xr2, Y1),
             i2c2=i2c2.os(X1, Y1) - i2c2.os(X2, Y1))

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
    # compute statistics of interest under the null
    D1.null <- discr.stat(Xn1, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr
    D2.null <- discr.stat(Xn2, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr

    icc1.null <- icc.os(Xnr1, Y1)
    icc2.null <- icc.os(Xnr2, Y1)
    i2c21.null <- i2c2.os(Xn1, Y1)
    i2c22.null <- i2c2.os(Xn2, Y1)
    return(list(discr=c(D1.null - D2.null, D2.null - D1.null),
                icc=c(icc1.null - icc2.null, icc2.null - icc1.null),
                i2c2=c(i2c21.null - i2c22.null, i2c22.null - i2c21.null)))
  }, mc.cores=no_cores)

  # compute null distribution of difference between discriminabilities
  null.diff <- list(discr=sapply(null.stats, function(x) x$discr),
                    icc=sapply(null.stats, function(x) x$icc),
                    i2c2=sapply(null.stats, function(x) x$i2c2))

  return(do.call(rbind, lapply(names(tr), function(stat.name) {
    data.frame(stat.name=stat.name, stat=tr[[stat.name]],
               p.value=mean(null.diff[[stat.name]] > tr[[stat.name]])*(nperm - 1)/nperm + 1/nperm)
    }))
  )
}

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
# a simulation where both pipelines are equally discriminable
# 2 classes
sim.no_signal <- function(n, d, sigma=0) {
  # samples 1 and 2 have classes are from same distribution, so no signal should be detected w.p. alpha
  samp1 <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  samp2 <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable, and pipeline 1 is more discriminable
# than pipeline 2
# 2 classes
sim.linear_sig <- function(n, d, sigma=0) {
  S <- diag(d)
  S[1, 1] <- 2
  S[-c(1), -c(1)] <- 1
  S2 <- S; S2[1,1] <- sigma  # sample 2 has greater covariance in signal dimension
  mus=cbind(rep(0, d), c(1, rep(0, d-1))) # with the same mean signal shift between the classes
  # sample 1 should be more discriminable than sample 2
  samp1 <- sim_gmm(mus=mus, Sigmas=abind(S, S, along=3), n)
  samp2 <- sim_gmm(mus=mus, Sigmas=abind(S, S, along=3), n)
  return(list(X1=samp1$X, X2=samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y))
}


## Crossed Signal Difference
# a simulation where classes are crossed but distinguishable
# also contains correlation btwn dimensions
# 2 classes
sim.crossed_sig <- function(n, d, K=16, sigma=0) {
  # class mus
  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)
  S.class <- diag(d)*sqrt(K)

  mus.class <- t(rbind(mvrnorm(n=K/2, mu.class.1, S.class),
                       mvrnorm(n=K/2, mu.class.2, S.class)))
  ni <- n/K

  # crossed signal
  Sigma.1 <- cbind(c(2,0), c(0,0.1))
  Sigma.2 <- cbind(c(0.1,0), c(0,2))
  mus=cbind(rep(0, d), rep(0, d))

  X1 <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    rho <- runif(1, min=-1, max=1)*sqrt(2*0.1)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rho
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rho
    sim <- sim_gmm(mus=mus, Sigmas=Sigmas, ni)
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  }))

  X2 <- do.call(rbind, lapply(1:K, function(k) {
    # add random correlation
    Sigmas <- abind(Sigma.1, Sigma.2, along = 3)
    rho <- runif(1, min=-1, max=1)*sqrt(2*0.1)
    Sigmas[1,2,1] <- Sigmas[2,1,1] <- rho
    Sigmas[1,2,2] <- Sigmas[2,1,2] <- -rho
    sim <- sim_gmm(mus=mus, Sigmas=Sigmas, ni)
    return(sweep(sim$X, 2, mus.class[,k], "+"))
  })) + array(rnorm(n*d), dim=c(n, d))*sigma

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X1=X1, X2=X2, Y=Y))
}

## Samples from Multiclass Gaussians
# a simulation where there are multiple classes present, and a correlation structure
# 2 classes
sim.multiclass_gaussian <- function(n, d, K=16, sigma=0) {
  S.k <- diag(d)*1
  S.k[upper.tri(S.k)] <- 0.5  # correlated
  S.k[lower.tri(S.k)] <- 0.5

  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K/2, mu.class.1, S.class),
                 mvrnorm(n=K/2, mu.class.2, S.class)))
  S1 <- abind(lapply(1:K, function(k) S.k), along=3)

  samp1 <- sim_gmm(mus=mus, Sigmas=S1, n)
  samp2 <- sim_gmm(mus=mus, Sigmas=S1, n)
  return(list(X1=samp1$X, X2= samp2$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp1$Y))
}

sim.multiclass_ann_disc <- function(n, d, K=16, sigma=0) {
  # centers
  K.cent <- K/2
  mu.class <- rep(0, d)
  S.class <- diag(d)*sqrt(K)

  mu.class.1 <- rep(0, d)
  mu.class.2 <- c(1, rep(0, d-1))*sqrt(K)
  S.class <- diag(d)*sqrt(K)

  mus <- t(rbind(mvrnorm(n=K.cent/2, mu.class.1, S.class),
                 mvrnorm(n=K.cent/2, mu.class.2, S.class)))

  ni <- n/K
  X1 <- do.call(rbind, lapply(1:K.cent, function(k) {
    X <- array(NaN, dim=c(ni*2, d))
    X[1:ni,] <- sweep(mgc.sims.2ball(ni, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni + 1):(2*ni),] <- sweep(mgc.sims.2sphere(ni, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  X2 <- do.call(rbind, lapply(1:K.cent, function(k) {
    X <- array(NaN, dim=c(ni*2, d))
    X[1:ni,] <- sweep(mgc.sims.2ball(ni, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni + 1):(2*ni),] <- sweep(mgc.sims.2sphere(ni, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X1=X1, X2=X2 + array(rnorm(n*d), dim=c(n, d))*sigma, Y=Y))
}

## --------------------------------------
# Driver
## --------------------------------------
n <- 128; d <- 2
nrep <- 500
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

list.results.ts <- mclapply(experiments, function(exper) {
  sim <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
  res <- test.two_sample(sim$X1, sim$X2, sim$Y)
  res$sim.name <- exper$sim.name; res$n <- n; res$d <- d; res$i <- exper$i
  res$sigma <- exper$sigma
  return(res)
}, mc.cores=no_cores)

ts.results <- do.call(rbind, list.results.ts)
saveRDS(list(ts.results=ts.results, list.results=list.results.ts), '../data/sims/discr_sims_ts.rds')

