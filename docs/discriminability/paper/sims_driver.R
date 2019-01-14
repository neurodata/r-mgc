# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
no_cores = detectCores() - 1

# redefine sothat the IO is the same for ICC/I2C2
iccadj <- function(Xr, Y) {
  if (!is.null(dim(Xr))) {
    if (dim(Xr)[2] > 1) {
      Xr <- lol.project.pca(Xr, r=1)$Xr
    }
  }
  data <- data.frame(x=Xr, y=Y)
  fit <- anova(aov(x ~ y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  a <- length(unique(Y))
  tmp.outj <- as.numeric(aggregate(x ~ y, data=data, FUN = length)$x)
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  p = fit[["Pr(>F)"]][1]
  return(list(srel=r, pval=p))
}

i2c2adj <- function(X, Y, Z=NULL) {
  if (is.null(Z)) {
    return(I2C2.original(X, Y, I=length(unique(Y)), twoway = FALSE)$lambda)
  } else {
    return(I2C2.original(X, Y, visit=Z)$lambda)
  }
}

manova <- function(X, Y) {
  fit <- manova(X ~ Y)
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  a <- length(unique(Y))
  tmp.outj <- as.numeric(aggregate(X ~ Y, data=data, FUN = length))
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  p = fit[["Pr(>F)"]][1]
}

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

nrep=50
n.max <- 512
n.min <- 16
d=2
rlen=10
opath='./data/sims'

log.seq <- function(from=0, to=30, length=rlen) {
  round(exp(seq(from=log(from), to=log(to), length.out=length)))
}

ns <- log.seq(from=n.min, to=n.max)

sims <- list(discr.sims.linear, discr.sims.linear, discr.sims.linear,# discr.sims.exp,
             discr.sims.cross, discr.sims.radial)#, discr.sims.beta)
sims.opts <- list(list(d=d, K=2, mean.scale=0, cov.scale=20), list(d=d, K=2, mean.scale=1, cov.scale=20),
                  list(d=d, K=5, mean.scale=1, cov.scale=20),#list(n=n, d=d, K=2, cov.scale=4),
                  list(d=d, K=2, cov.scale=20), list(d=d, K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", "Linear, 5 Class", "Cross", "Radial")

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
  mag <- apply(X, c(1), dist)
  ang <- apply(X, c(1), function(x) {
    uvec <- array(0, length(x)); uvec[1] <- 1
    return(acos(x %*% uvec/(sqrt(sum(x^2)))))})
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
