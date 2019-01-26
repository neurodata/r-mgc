# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
no_cores = detectCores() - 1

# -----------------------------
# Simulations
# -----------------------------

# redefine simulations so that we can obtain the best/worst dimensions relatively
# easily
sim.linear.os <- function(...) {
  simout <- discr.sims.linear(...)
  simout$d.best <- simout$X[,1]  # first dimension has all signal
  simout$d.worst <- simout$X[,2]  # second has none
  return(simout)
}

sim.cross.os <- function(...) {
  simout <- discr.sims.cross(...)
  simout$d.best <- simout$X[,1]  # either first or second dimension are top signal wise
  simout$d.worst <- simout$X %*% array(c(1, 1), dim=c(dim(simout$X)[2])) # worst dimension is the direction of maximal variance
  return(simout)
}

sim.radial.os <- function(...) {
  simout <- discr.sims.radial(...)
  simout$d.best <- apply(simout$X, c(1), dist)  # best dimension is the radius of the point
  # worst dimension is the angle of the point in radians relative 0 rad
  simout$d.worst <- apply(simout$X, c(1), function(x) {
    uvec <- array(0, length(x)); uvec[1] <- 1
    return(acos(x %*% uvec/(sqrt(sum(x^2)))))})
  return(simout)
}

# simulations and options as a list
sims <- list(sim.linear.os, sim.linear.os, #sim.linear,# discr.sims.exp,
             sim.cross.os, sim.radial.os)#, discr.sims.beta)
sims.opts <- list(list(K=2, signal.lshift=0), list(K=2, signal.lshift=1),
                  #list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(K=2, cov.scale=20), list(K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", #"Linear, 5 Class",
                "Cross", "Radial")

# -----------------------------
# One-Sample Testing Algorithms
# -----------------------------
# define all the relevant one-sample tests with a similar interface for simplicity in the driver

# one-way anova
anova.os <- function(x, y) {
  data <- data.frame(x=x, y=y)
  fit <- anova(aov(x ~ y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  f = fit[["F value"]][1]
  p = fit[["Pr(>F)"]][1]
  return(list(f=f, p=p))
}

# driver for one-way anova
anova.onesample.driver <- function(sim, ...) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("ANOVA, best", "ANOVA, worst")
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    res <- anova.os(Xrs[[i]], sim$Y)
    return(data.frame(alg=algs[[i]], tstat=res$f, pval=res$p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# one-way ICC
icc.os <- function(x, y) {
  data <- data.frame(x=x, y=y)
  fit <- anova(aov(x ~ y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  a <- length(unique(sim$Y))
  tmp.outj <- as.numeric(aggregate(x ~ y, data=data, FUN = length)$x)
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  return(r)
}

# driver for one-way ICC
icc.onesample.driver <- function(sim, nperm=100, ...) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("ICC, best", "ICC, worst")
  N <- length(sim$Y)
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    # relative statistic
    r <- icc.os(Xrs[[i]], sim$Y)
    # permutation approach for p-value
    nr <- sapply(1:nperm, function(j) {
      samplen <- sample(N)
      # randomly permute labels
      perm.stat <- do.call(icc.os, list(x=Xrs[[i]], y=sim$Y[samplen]))
      return(perm.stat)
    })
    # p-value is fraction of times statistic of permutations
    # is more extreme than the relative statistic
    p <- (sum(nr>r) + 1)/(nperm + 1)
    return(data.frame(alg=algs[[i]], tstat=r, pval=p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# one-sample MANOVA
manova.onesample.driver <- function(sim) {
  fit <- manova(sim$X ~ sim$Y)
  return(data.frame(alg="MANOVA", tstat=summary(fit)$stats["sim$Y", "approx F"],
                    pval=summary(fit)$stats["sim$Y", "Pr(>F)"]))
}

# I2C2 wrapper
i2c2.os <- function(X, Y) {
  return(I2C2.original(y=X, id=Y, visit=rep(1, length(Y)), twoway=FALSE)$lambda)
}

# one-sample I2C2
i2c2.onesample.driver <- function(sim, nperm=100, ...) {
  # relative statistic
  N <- length(sim$Y)
  r <- i2c2.os(sim$X, sim$Y)
  # permutation approach for p-value
  nr <- sapply(1:nperm, function(i) {
    samplen <- sample(N)
    # randomly permute labels
    do.call(i2c2.os, list(X=sim$X, Y=sim$Y[samplen]))
  })
  # p-value is fraction of times statistic of permutations
  # is more extreme than the relative statistic
  p <- (sum(nr>r) + 1)/(nperm + 1)
  return(data.frame(alg="I2C2", tstat=r, pval=p))
}

discr.onesample.driver <- function(sim, nperm=100, ...) {
  # relative statistic
  N <- length(sim$Y)
  D <- discr.distance(sim$X)
  r <- discr.stat(D, sim$Y, is.dist=TRUE)$discr
  # permutation approach for p-value
  nr <- sapply(1:nperm, function(i) {
    samplen <- sample(N)
    # randomly permute labels
    do.call(discr.stat, list(X=D, Y=sim$Y[samplen], is.dist=TRUE))$discr
  })
  # p-value is fraction of times statistic of permutations
  # is more extreme than the relative statistic
  p <- (sum(nr>r) + 1)/(nperm + 1)
  return(data.frame(alg="Discr", tstat=r, pval=p))
}

algs <- list(discr.onesample.driver, anova.onesample.driver, icc.onesample.driver,
             manova.onesample.driver, i2c2.onesample.driver)
# ----------------------------------
## One-Sample Driver
# ----------------------------------

# Set Global Parameters for Investigating
nrep=200  # number of iterations per n, trial
n.max <- 512  # maximum number of samples
n.min <- 16  # minimum number of samples
d=2  # number of dimensions
rlen=10  # number of ns to try
opath='./data/sims'  # output path

log.seq <- function(from=0, to=30, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

ns <- log.seq(from=n.min, to=n.max)

experiments <- do.call(c, lapply(seq_along(sims), function(sims, sims.names, sims.opts, i) {
  sim <- sims[[i]]; sim.name <- sims.names[[i]]; sim.opts <- sims.opts[[i]]
  do.call(c, lapply(ns, function(n) {
    do.call(c, lapply(1:nrep, function(j) {
      sim.out <- do.call(sim, c(sim.opts, n=n, d=d))
      lapply(algs, function(alg) {
        return(list(sim.name=sim.name, sim=sim.out, i=j, n=n, d=d, alg=alg))
      })
    }))
  }))
}, sims=sims, sims.names=sims.names, sims.opts=sims.opts))

## One Sample Results
# mcapply over the number of repetitions
results <- mclapply(experiments, function(exp) {
  res <- do.call(exp$alg, list(exp$sim))
  return(data.frame(sim.name=exp$sim.name, n=exp$n, i=exp$i, alg=res$alg, tstat=res$tstat, pval=res$pval))
}, mc.cores=no_cores)

results <- do.call(rbind, results)
saveRDS(results, file.path(opath, paste('discr_sims_os', '.rds', sep="")))

# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
no_cores = detectCores() - 1

# -----------------------------
# TS Simulations
# -----------------------------

# redefine simulations so that we can obtain the best/worst dimensions relatively
# easily
sim.linear.ts <- function(...) {
  simout <- discr.sims.linear(...)
  X.g2 <- simout$X + array(rnorm(prod(dim(simout$X))), dim=dim(simout$X))
  simout$X <- rbind(simout$X, X.g2)
  simout$Z <- factor(c(rep(1, length(simout$Y)), rep(2, length(simout$Y))))
  simout$Y <- c(simout$Y, simout$Y)
  simout$d.best <- simout$X[,1]  # first dimension has all signal
  simout$d.worst <- simout$X[,2]  # second has none
  return(simout)
}

sim.cross.ts <- function(...) {
  simout <- discr.sims.cross(...)
  X.g2 <- simout$X + array(rnorm(prod(dim(simout$X))), dim=dim(simout$X))
  simout$X <- rbind(simout$X, X.g2)
  simout$Z <- factor(c(rep(1, length(simout$Y)), rep(2, length(simout$Y))))
  simout$Y <- c(simout$Y, simout$Y)
  simout$d.best <- simout$X[,1]  # either first or second dimension are top signal wise
  simout$d.worst <- simout$X %*% array(c(1, 1), dim=c(dim(simout$X)[2])) # worst dimension is the direction of maximal variance
  return(simout)
}

sim.radial.ts <- function(...) {
  simout <- discr.sims.radial(...)
  X.g2 <- simout$X + array(rnorm(prod(dim(simout$X))), dim=dim(simout$X))
  simout$X <- rbind(simout$X, X.g2)
  simout$Z <- factor(c(rep(1, length(simout$Y)), rep(2, length(simout$Y))))
  simout$Y <- c(simout$Y, simout$Y)
  simout$d.best <- apply(simout$X, c(1), dist)  # best dimension is the radius of the point
  # worst dimension is the angle of the point in radians relative 0 rad
  simout$d.worst <- apply(simout$X, c(1), function(x) {
    uvec <- array(0, length(x)); uvec[1] <- 1
    return(acos(x %*% uvec/(sqrt(sum(x^2)))))})
  return(simout)
}

# simulations and options as a list
sims <- list(sim.linear.ts, sim.linear.ts, #sim.linear,# discr.sims.exp,
             sim.cross.ts, sim.radial.ts)#, discr.sims.beta)
sims.opts <- list(list(K=2, signal.lshift=0), list(K=2, signal.lshift=1),
                  #list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(K=2, signal.scale=20), list(K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", #"Linear, 5 Class",
                "Cross", "Radial")

# -----------------------------
# Two-Sample Testing Algorithms
# -----------------------------
# define all the relevant two-sample tests with a similar interface for simplicity in the driver

# two-way anova
anova.ts <- function(x, y, z) {
  data <- data.frame(x=x, y=y, z=z)
  fit <- anova(aov(x ~ y*z, data=data))
  f = fit[["F value"]][3]
  p = fit[["Pr(>F)"]][3]
  return(list(f=f, p=p))
}

# driver for two-way anova
anova.twosample.driver <- function(sim, ...) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("ANOVA, best", "ANOVA, worst")
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    res <- anova.ts(Xrs[[i]], sim$Y, sim$Z)
    return(data.frame(alg=algs[[i]], tstat=res$f, pval=res$p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# two-way ICC
icc.ts <- function(x, y, z) {
  data <- data.frame(x=x, y=y, z=z)
  fit <- anova(aov(x ~ y*z, data=data))
  MSa <- fit$"Mean Sq"[3]
  MSw <- var.w <- fit$"Mean Sq"[3]
  a <- length(unique(sim$Y))
  tmp.outj <- as.numeric(aggregate(x ~ y, data=data, FUN = length)$x)
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  return(r)
}

# driver for two-way ICC
icc.twosample.driver <- function(sim, nperm=100, ...) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("ICC, best", "ICC, worst")
  N <- length(sim$Y)
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    # relative statistic
    r <- icc.os(Xrs[[i]], sim$Y)
    # permutation approach for p-value
    nr <- sapply(1:nperm, function(j) {
      samplen <- sample(N)
      # randomly permute labels
      perm.stat <- do.call(icc.ts, list(x=Xrs[[i]], y=sim$Y[samplen], z=sim$Z[samplen]))
      return(perm.stat)
    })
    # p-value is fraction of times statistic of permutations
    # is more extreme than the relative statistic
    p <- (sum(nr>r) + 1)/(nperm + 1)
    return(data.frame(alg=algs[[i]], tstat=r, pval=p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# two-sample MANOVA
manova.twosample.driver <- function(sim) {
  fit <- manova(sim$X ~ sim$Y*sim$Z)
  return(data.frame(alg="MANOVA", tstat=summary(fit)$stats["sim$Y:sim$Z", "approx F"],
                    pval=summary(fit)$stats["sim$Y:sim$Z", "Pr(>F)"]))
}

# I2C2 wrapper
i2c2.ts <- function(X, Y, Z) {
  return(I2C2.original(y=X, id=Y, visit=Z, twoway=TRUE)$lambda)
}

# two-way I2C2
i2c2.twosample.driver <- function(sim, nperm=100, ...) {
  # relative statistic
  N <- length(sim$Y)
  r <- i2c2.ts(sim$X, sim$Y, sim$Z)
  # permutation approach for p-value
  nr <- sapply(1:nperm, function(i) {
    samplen <- sample(N)
    # randomly permute labels
    do.call(i2c2.ts, list(X=sim$X, Y=sim$Y[samplen], Z=sim$Z[samplen]))
  })
  # p-value is fraction of times statistic of permutations
  # is more extreme than the relative statistic
  p <- (sum(nr>r) + 1)/(nperm + 1)
  return(data.frame(alg="I2C2", tstat=r, pval=p))
}

discr.twosample.driver <- function(sim, nperm=100, ...) {
  D1 <- discr.distance(sim$X[sim$Z == 1,])
  D2 <- discr.distance(sim$X[sim$Z == 2,])
  res <- discr.test.two_sample(D1, D2, ids=sim$Y[sim$Z == 1])
  return(data.frame(alg="Discr", pval=res$pval))
}

algs <- list(discr.twosample.driver, anova.twosample.driver, #icc.twosample.driver,
             manova.twosample.driver)#, i2c2.twosample.driver)
# ----------------------------------
## Two-Sample Driver
# ----------------------------------

# Set Global Parameters for Investigating
nrep=50  # number of iterations per n, trial
n.max <- 512  # maximum number of samples
n.min <- 16  # minimum number of samples
d=2  # number of dimensions
rlen=10  # number of ns to try
opath='./data/sims'  # output path

log.seq <- function(from=0, to=30, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

ns <- log.seq(from=n.min, to=n.max)

experiments <- do.call(c, lapply(seq_along(sims), function(sims, sims.names, sims.opts, i) {
  sim <- sims[[i]]; sim.name <- sims.names[[i]]; sim.opts <- sims.opts[[i]]
  do.call(c, lapply(ns, function(n) {
    do.call(c, lapply(1:nrep, function(j) {
      sim.out <- do.call(sim, c(sim.opts, n=n, d=d))
      lapply(algs, function(alg) {
        return(list(sim.name=sim.name, sim=sim.out, i=j, n=n, d=d, alg=alg))
      })
    }))
  }))
}, sims=sims, sims.names=sims.names, sims.opts=sims.opts))

## Two Sample Results
# mcapply over the number of repetitions
results <- mclapply(experiments, function(exp) {
  res <- do.call(exp$alg, list(exp$sim))
  return(data.frame(sim.name=exp$sim.name, n=exp$n, i=exp$i, alg=res$alg, pval=res$pval))
}, mc.cores=no_cores)

results <- do.call(rbind, results)
saveRDS(results, file.path(opath, paste('discr_sims_ts', '.rds', sep="")))
