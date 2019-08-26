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
                  list(K=2, signal.scale=5, mean.scale=1), list(K=2))#,
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
    return(data.frame(alg=algs[[i]], tstat=res$f, p.value=res$p))
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
  a <- length(unique(y))
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
    return(data.frame(alg=algs[[i]], tstat=r, p.value=p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# one-sample MANOVA
manova.onesample.driver <- function(sim) {
  fit <- manova(sim$X ~ sim$Y)
  return(data.frame(alg="MANOVA", tstat=summary(fit)$stats["sim$Y", "approx F"],
                    p.value=summary(fit)$stats["sim$Y", "Pr(>F)"]))
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
  return(data.frame(alg="I2C2", tstat=r, p.value=p))
}

discr.onesample.driver <- function(sim, nperm=100, ...) {
  res=discr.test.one_sample(sim$X, sim$Y)
  return(data.frame(alg="Discr", tstat=res$stat, p.value=res$p.value))
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
  tryCatch({
    res <- do.call(exp$alg, list(exp$sim))
    return(data.frame(sim.name=exp$sim.name, n=exp$n, i=exp$i, alg=res$alg, tstat=res$tstat, p.value=res$p.value))
  }, error=function(e) {return(NULL)})
}, mc.cores=no_cores)

results <- do.call(rbind, results)
saveRDS(results, file.path(opath, paste('discr_sims_os', '.rds', sep="")))

# Parallelize Stuff
#=========================#
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
sim.linear.ts <- function(opt1, opt2, n, d) {
  s.g1 <- do.call(discr.sims.linear, c(opt1, list(n=n, d=d)))
  s.g2 <- do.call(discr.sims.linear, c(opt2, list(n=n*3, d=d)))
  g2.out <- list(X=NULL, Y=NULL)
  for (y in unique(s.g1$Y)) {
    idx.g2 <- which(s.g2$Y == y)
    n.y <- sum(s.g1$Y == y)
    g2.out$X <- rbind(g2.out$X, s.g2$X[idx.g2[1:n.y],])
    g2.out$Y <- c(g2.out$Y, s.g2$Y[idx.g2[1:n.y]])
  }
  ord.g1 <- order(s.g1$Y)
  simout <- list(X=rbind(s.g1$X[ord.g1,], g2.out$X),
                 Z=factor(c(rep(1, n), rep(2, n))),
                 Y=c(s.g1$Y[ord.g1], g2.out$Y))
  simout$d.best <- simout$X[,1]  # first dimension has all signal
  simout$d.worst <- simout$X[,2]  # second has none
  return(simout)
}

sim.cross.ts <- function(opt1, opt2, n, d) {
  s.g1 <- do.call(discr.sims.cross, c(opt1, list(n=n, d=d)))
  s.g2 <- do.call(discr.sims.cross, c(opt2, list(n=n*3, d=d)))
  g2.out <- list(X=NULL, Y=NULL)
  for (y in unique(s.g1$Y)) {
    idx.g2 <- which(s.g2$Y == y)
    n.y <- sum(s.g1$Y == y)
    g2.out$X <- rbind(g2.out$X, s.g2$X[idx.g2[1:n.y],])
    g2.out$Y <- c(g2.out$Y, s.g2$Y[idx.g2[1:n.y]])
  }
  ord.g1 <- order(s.g1$Y)
  simout <- list(X=rbind(s.g1$X[ord.g1,], g2.out$X),
                 Z=factor(c(rep(1, n), rep(2, n))),
                 Y=c(s.g1$Y[ord.g1], g2.out$Y))
  simout$d.best <- simout$X[,1]  # either first or second dimension are top signal wise
  simout$d.worst <- simout$X %*% array(c(1, 1), dim=c(dim(simout$X)[2])) # worst dimension is the direction of maximal variance which is the diagonal
  return(simout)
}

sim.radial.ts <- function(opt1, opt2, n, d) {
  s.g1 <- do.call(discr.sims.radial, c(opt1, list(n=n, d=d)))
  s.g2 <- do.call(discr.sims.radial, c(opt2, list(n=n*3, d=d)))
  g2.out <- list(X=NULL, Y=NULL)
  for (y in unique(s.g1$Y)) {
    idx.g2 <- which(s.g2$Y == y)
    n.y <- sum(s.g1$Y == y)
    g2.out$X <- rbind(g2.out$X, s.g2$X[idx.g2[1:n.y],])
    g2.out$Y <- c(g2.out$Y, s.g2$Y[idx.g2[1:n.y]])
  }
  ord.g1 <- order(s.g1$Y)
  simout <- list(X=rbind(s.g1$X[ord.g1,], g2.out$X),
                 Z=factor(c(rep(1, n), rep(2, n))),
                 Y=c(s.g1$Y[ord.g1], g2.out$Y))
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
sims.opts <- list(list(opt1=list(K=2, signal.lshift=0), opt2=list(K=2, signal.lshift=0)),
                  list(opt1=list(K=2, signal.lshift=1), opt2=list(K=2, signal.lshift=1, signal.scale=2)),
                  #list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(opt1=list(K=2, signal.scale=5, mean.scale=1), opt2=list(K=2, signal.scale=5, non.scale=5, mean.scale=1)),
                  list(opt1=list(K=2), opt2=list(K=2, er.scale=0.2)))#,
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
  f = fit["y:z", "F value"]
  p = fit["y:z", "Pr(>F)"]
  return(list(f=f, p=p))
}

# driver for two-way anova
anova.twosample.driver <- function(sim, ...) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("ANOVA, best", "ANOVA, worst")
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    res <- anova.ts(Xrs[[i]], sim$Y, sim$Z)
    return(data.frame(alg=algs[[i]], tstat=res$f, p.value=res$p))
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
    return(data.frame(alg=algs[[i]], tstat=r, p.value=p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

# two-sample MANOVA
manova.twosample.driver <- function(sim) {
  fit <- manova(sim$X ~ sim$Y*sim$Z)
  return(data.frame(alg="MANOVA", tstat=summary(fit)$stats["sim$Y:sim$Z", "approx F"],
                    p.value=summary(fit)$stats["sim$Y:sim$Z", "Pr(>F)"]))
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
  return(data.frame(alg="I2C2", tstat=r, p.value=p))
}

discr.twosample.driver <- function(sim, nperm=100, ...) {
  X1 <- discr.distance(sim$X[sim$Z == 1,])
  X2 <- discr.distance(sim$X[sim$Z == 2,])
  res <- discr.test.two_sample(X1, X2, Y=sim$Y[sim$Z == 1])
  return(data.frame(alg="Discr", p.value=res$p.value))
}

algs <- list(discr.twosample.driver, anova.twosample.driver, #icc.twosample.driver,
             manova.twosample.driver)#, i2c2.twosample.driver)
# ----------------------------------
## Two-Sample Driver
# ----------------------------------

# Set Global Parameters for Investigating
nrep=200  # number of iterations per n, trial
n.max <- 512  # maximum number of samples
n.min <- 16  # minimum number of samples
d=2  # number of dimensions
rlen=10  # number of ns to try
opath='./data/sims'  # output path
ts.sim.opath='./data/sims/ts_sims'

log.seq <- function(from=0, to=30, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

ns <- log.seq(from=n.min, to=n.max)

experiments <- do.call(c, lapply(seq_along(sims), function(sims, sims.names, sims.opts, i) {
  sim <- sims[[i]]; sim.name <- sims.names[[i]]; sim.opts <- sims.opts[[i]]
  do.call(c, lapply(ns, function(n) {
    do.call(c, lapply(1:nrep, function(j) {
      sim.out <- do.call(sim, list(opt1=sim.opts$opt1, opt2=sim.opts$opt2, n=n, d=d))
      lapply(algs, function(alg) {
        return(list(sim.name=sim.name, sim=sim.out, i=j, n=n, d=d, alg=alg))
      })
    }))
  }))
}, sims=sims, sims.names=sims.names, sims.opts=sims.opts))

dir.create(ts.sim.opath)
## Two Sample Results
# mcapply over the number of repetitions
results <- mclapply(experiments, function(exp) {
 tryCatch({
    res <- do.call(exp$alg, list(exp$sim))
    return(data.frame(sim.name=exp$sim.name, n=exp$n, i=exp$i, alg=res$alg, p.value=res$p.value))
  }, error=function(e){return(NULL)})
}, mc.cores=no_cores)

results <- do.call(rbind, results)
saveRDS(results, file.path(opath, paste('discr_sims_ts', '.rds', sep="")))


