##----------------------------------------------
# One Sample Tests
##----------------------------------------------
require(lolR)
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


## One Sample Driver
test.one_sample <- function(X, Y, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                            dist.return=NULL, remove.isolates=TRUE, nperm=100, no_cores=1) {

  validated <- mgc:::discr.validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                     remove.isolates=remove.isolates)

  D <- validated$D; Y <- validated$Y; N <- nrow(D)
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  } else if (no_cores >= 0.8*detectCores()) {
    warning("You have requested a number of cores near your machine's core count. Expected decreased performance.")
  }
  Xr <- lol.project.pca(X, 1)$Xr
  # compute references for the statistics
  tr <- list(discr=discr.stat(D, Y, is.dist=TRUE)$discr,
             icc=icc.os(Xr, Y),
             i2c2=i2c2.os(X, Y))

  nr <- mclapply(1:nperm, function(i) {
    perm.Y <- Y[sample(N)]
    return(list(discr=discr.stat(D, perm.Y, is.dist=TRUE)$discr,
                icc=icc.os(Xr, perm.Y),
                i2c2=i2c2.os(X, perm.Y)))
  }, mc.cores=no_cores)

  null.stats <- list(discr=lapply(nr, function(x) x$discr),
                     icc=lapply(nr, function(x) x$icc),
                     i2c2=lapply(nr, function(x) x$i2c2))

  return(do.call(rbind, lapply(names(tr), function(stat.name) {
    data.frame(stat.name=stat.name, stat=tr[[stat.name]],
               p.value=mean(null.stats[[stat.name]] > tr[[stat.name]])*(nperm - 1)/nperm + 1/nperm)
  })))
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
# a simulation where no distinguishable signal present
# 2 classes
sim.no_signal <- function(n, d, sigma=1) {
  # classes are from same distribution, so signal should be detected w.p. alpha
  samp <- sim_gmm(mus=cbind(rep(0, d), rep(0,d)), Sigmas=abind(diag(d), diag(d), along=3), n)
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y))
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
  mus=t(mvrnorm(n=2, c(0, 0), S.class)) # with a mean signal shift between the classes
  samp <- sim_gmm(mus=mus, Sigmas=abind(Sigma, Sigma, along=3), n)
  return(list(X=samp$X + array(rnorm(n*d), dim=c(n, d))*sigma, Y=samp$Y))
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
  return(list(X=X, Y=Y))
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
  Sigmas <- abind(lapply(1:K, function(k) S.k), along=3)

  samp <- sim_gmm(mus=mus, Sigmas=Sigmas, n)
  return(list(X=samp$X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=samp$Y))
}

# 8 pairs of annulus/discs
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

  X <- do.call(rbind, lapply(1:K.cent, function(k) {
    X <- array(NaN, dim=c(ni*2, d))
    X[1:ni,] <- sweep(mgc.sims.2ball(ni, d, r=1, cov.scale=0.1), 2, mus[,k], "+")
    X[(ni + 1):(2*ni),] <- sweep(mgc.sims.2sphere(ni, r=1, d=d, cov.scale=0.1), 2, mus[,k], "+")
    return(X)
  }))

  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni)))
  return(list(X=X + array(rnorm(n*d)*sigma, dim=c(n, d)), Y=Y))
}

## -------------------------
# Driver
## -------------------------
n <- 128; d <- 2
nrep <- 500
n.sigma <- 15

simulations <- list(sim.no_signal, sim.linear_sig, sim.crossed_sig,
                    sim.multiclass_gaussian, sim.multiclass_ann_disc)
sims.sig.max <- c(20, 20, 20, 20, 20)
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


list.results.os <- mclapply(experiments, function(exper) {
  sim <- do.call(exper$sim, list(n=n, d=d, sigma=exper$sigma))
  res <- test.one_sample(sim$X, sim$Y)
  res$sim.name <- exper$sim.name; res$n <- n; res$d <- d; res$i <- exper$i
  res$sigma <- exper$sigma
  return(res)
}, mc.cores=no_cores)

bound.results.os <- do.call(rbind, list.results.os)
saveRDS(list(os.results=bound.results.os, list.results=list.results.os), '../data/sims/discr_sims_os.rds')

