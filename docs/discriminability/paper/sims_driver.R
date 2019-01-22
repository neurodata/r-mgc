# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
no_cores = detectCores() - 1

# Set Global Parameters for Investigating
nrep=50  # number of iterations per n, trial
n.max <- 512  # maximum number of samples
n.min <- 16  # minimum number of samples
d=2  # number of dimensions
rlen=10  # number of ns to try
opath='./data/sims'  # output path


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
  simout$d.best <- apply(X, c(1), dist)  # best dimension is the radius of the point
  # worst dimension is the angle of the point in radians relative 0 rad
  simout$d.worst <- apply(X, c(1), function(x) {
    uvec <- array(0, length(x)); uvec[1] <- 1
    return(acos(x %*% uvec/(sqrt(sum(x^2)))))})
}

sims <- list(sim.linear.os, sim.linear.os, #sim.linear,# discr.sims.exp,
             sim.cross.os, sim.radial.os)#, discr.sims.beta)
sims.opts <- list(list(d=d, K=2, signal.lshift=0), list(d=d, K=2, signal.lshift=1),
                  #list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(d=d, K=2, cov.scale=20), list(d=d, K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", #"Linear, 5 Class",
                "Cross", "Radial")

anova.onesample <- function(sim) {
  Xrs <- list(best=sim$d.best, worst=sim$d.worst)
  names(Xrs) <- c("Anova, best", "Anova, worst")
  result <- lapply(seq_along(Xrs), function(Xrs, algs, i) {
    data <- data.frame(x=Xrs[[i]], y=sim$Y)
    fit <- anova(aov(x ~ y, data=data))
    MSa <- fit$"Mean Sq"[1]
    MSw <- var.w <- fit$"Mean Sq"[2]
    a <- length(unique(sim$Y))
    tmp.outj <- as.numeric(aggregate(x ~ y, data=data, FUN = length)$x)
    k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
    var.a <- (MSa - MSw)/k
    r <- var.a/(var.w + var.a)
    p = fit[["Pr(>F)"]][1]
    return(data.frame(alg=algs[[i]], icc=r, tstat=r, pval=p))
  }, Xrs=Xrs, algs=names(Xrs))
  res <- do.call(rbind, result)
  return(res)
}

manova.onesample <- function(sim) {
  fit <- manova(sim$X ~ sim$Y)
  return(data.frame("Manova", i2c2=, i2c2.pval=, tstat=summary(fit)$stats["sim$Y", "approx F"],
                    pval=summary(fit)$stats["sim$Y", "Pr(>F)"]))
}

discr.onesample <- function()

one.sample.i2c2 <- function(X, ids, Z=NULL, metr=i2c2adj, nperm=100, verbose=FALSE) {
  N <- length(ids)
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  tr <- do.call(i2c2adj, list(X, ids))

  nr <- rep(0,nperm)
  for (i in 1:nperm){
    if (verbose) {
      print(i)
    }
    samplen <- sample(N)
    if (is.null(Z)) {
      Z = Z[samplen]
    }
    nr[i] <- do.call(i2c2adj, list(X=X[samplen,], Y=ids[samplen], Z=Z))
  }
  result <- list()
  result$srel <- tr
  result$null <- sort(nr)
  result$pval <- (sum(nr>tr) + 1)/(nperm + 1)
  return(result)
}

log.seq <- function(from=0, to=30, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

ns <- log.seq(from=n.min, to=n.max)

one.sample.tests <- list(discr.test.one_sample, iccadj, one.sample.i2c2)
stat.names <- c("Discr", "ICC", "I2C2")

## One Sample Results
# mcapply over the number of repetitions
results <- mclapply(1:nrep, function(i) {
  results <- data.frame(sim=c(), iteration=c(), stat=c(), n=c(), tstat=c(), pval=c())
  for (j in 1:length(sims)) {
    print(sprintf("j: %d", j))
    sim <- sims[[j]]; sim.opt <- sims.opts[[j]]; sim.name <- sims.names[[j]]
    for (k in 1:length(ns)) {
      n <- ns[k]
      print(sprintf("n: %d", n))
      sim.opt$n <- n
      sim.dat <- do.call(sim, sim.opt); X <- sim.dat$X; Y <- sim.dat$Y
      for (l in 1:length(one.sample.tests)) {
        ost <- one.sample.tests[[l]]; stat.name <- stat.names[[l]]
        t.res <- do.call(ost, list(X, Y)); tstat <- t.res$srel; pval <- t.res$pval
        results <- rbind(results,  data.frame(sim=sim.name, iteration=i, stat=stat.name, n=n, tstat=tstat, pval=pval))
      }
    }
  }
  return(results)
}, mc.cores=no_cores)

results <- do.call(rbind, results)
saveRDS(results, file.path(opath, paste('discr_sims_os', '.rds', sep="")))

two.sample.pca <- function(X, Y) {
  pca.res <- lol.project.pca(X, r=2)
  Xr <- pca.res$Xr
  return(list(Y=Y, X1=Xr[, 1], X2=Xr[, 2]))
}

two.sample.spherical <- function(X, Y) {

  return(list(Y=Y, X1=mag, X2=ang))
}

two.sample.tests <- list(discr.test.two_sample, iccadj)
stat.names <- c("Discr", "ICC")

sims <- list(discr.sims.linear, discr.sims.linear, discr.sims.linear,# discr.sims.exp,
             discr.sims.cross, discr.sims.radial)#, discr.sims.beta)
sims.opts <- list(list(d=d, K=2, mean.scale=0, cov.scale=20), list(d=d, K=2, mean.scale=1, cov.scale=20),
                  list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(d=d, K=2, cov.scale=20), list(d=d, K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", "Linear, 5 Class", "Cross", "Radial")
sim.twos_prep <- list(two.sample.pca, two.sample.pca, two.sample.pca, two.sample.pca, two.sample.spherical)
## Two Sample Results
# mcapply over the number of repetitions
results <- mclapply(1:nrep, function(i) {
  results <- data.frame(sim=c(), iteration=c(), stat=c(), n=c(),tstat=c(), pval=c())
  for (j in 1:length(sims)) {
    print(sprintf("j: %d", j))
    sim <- sims[[j]]; sim.opt <- sims.opts[[j]]; sim.name <- sims.names[[j]]
    for (k in 1:length(ns)) {
      n <- ns[k]
      print(sprintf("n: %d", n))
      sim.opt$n <- n
      sim.dat <- do.call(sim, sim.opt); X <- sim.dat$X; Y <- sim.dat$Y
      twosample.dat <- do.call(sim.twos_prep[[j]], list(X, Y))
      X1 <- twosample.dat$X1; X2 <- twosample.dat$X2; Xcomb <- c(twosample.dat$X1, twosample.dat$X2); Y <- Y; Ycomb <- c(Y, Y); Zcomb <- c(rep(1, length(Y)), rep(2, length(Y)))
      for (l in 1:length(two.sample.tests)) {
        tryCatch({
          tst <- two.sample.tests[[l]]; stat.name <- stat.names[[l]]
          if (stat.name == "Discr") {
            D1 <- discr.distance(X1); D2 <- discr.distance(X2)
            t.res <- do.call(tst, list(D1, D2, Y)); pval <- t.res$pval
          } else if (stat.name == "ICC") {
            t.res <- do.call(tst, list(Xcomb, as.factor(Ycomb))); pval <- t.res$pval
          } else {
            t.res <- do.call(tst, list(X=Xcomb, Y=as.factor(Ycomb), Z=Zcomb)); pval <- t.res$pval
          }
          results <- rbind(results,  data.frame(sim=sim.name, iteration=i, stat=stat.name, n=n, pval=pval))},
          error=function(e) print(e))
      }
    }
  }
  return(results)
}, mc.cores=no_cores)

results <- do.call(rbind, results)

saveRDS(results, file.path(opath, paste('discr_sims_ts', '.rds', sep="")))
