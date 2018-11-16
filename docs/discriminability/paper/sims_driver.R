# Parallelize Stuff
#=========================#
require(lolR)
require(MASS)
library(parallel)
no_cores = detectCores() - 1

# redefine sothat the IO is the same for ICC/I2C2
iccadj <- function(X, Y) {
  return(ICCbare(Y, lol.project.pca(X, r=1)$Xr))
}

i2c2adj <- function(X, Y) {
  return(i2c2(X, Y, visit=rep(1, length(Y)))$lambda)
}

one.sample.icc_i2c2 <- function(X, ids, metric, nperm=100, verbose=FALSE) {
  N <- length(ids)
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  tr <- do.call(metric, list(X, ids))

  nr <- rep(0,nperm)
  for (i in 1:nperm){
    if (verbose) {
      print(i)
    }
    samplen <- sample(N)
    nr[i] <- do.call(metric, list(X[samplen,], ids[samplen]))
  }
  result <- list()
  result$srel <- tr
  result$null <- sort(nr)
  result$pval <- (sum(nr>tr) + 1)/(nperm + 1)
  return(result)
}

one.sample.icc <- function(X, ids, nperm=100, verbose=FALSE) {
  one.sample.icc_i2c2(X, ids, iccadj, nperm=nperm, verbose=verbose)
}

one.sample.i2c2 <- function(X, ids, nperm=100, verbose=FALSE) {
  one.sample.icc_i2c2(X, ids, i2c2adj, nperm=nperm, verbose=verbose)
}


nrep=20
n.max <- 1024
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
                  list(d=d, K=2, cov.scale=5), list(d=d, K=2))#,
sims.names <- c("No Signal", "Linear, 2 Class", "Linear, 5 Class", "Cross", "Radial")

one.sample.tests <- list(discr.test.one_sample, one.sample.icc, one.sample.i2c2)
stat.names <- c("Discr", "ICC", "I2C2")

## One Sample Results
# mcapply over the number of repetitions
results <- mclapply(1:nrep, function(i) {
  results <- data.frame(sim=c(), iteration=c(), stat=c(), n=c(), pval=c())
  for (j in length(sims)) {
   sim <- sims[[j]]; sim.opt <- sims.opts[[j]]; sim.name <- sims.names[[j]]
    for (k in 1:length(ns)) {
      n <- ns[k]
      sim.opt$n <- n
      sim.dat <- do.call(sim, sim.opt); X <- sim.dat$X; Y <- sim.dat$Y
      for (l in 1:length(one.sample.tests)) {
        ost <- one.sample.tests[[l]]; stat.name <- stat.names[[l]]
        t.res <- do.call(ost, list(X, Y)); tstat <- t.res$srel; pval <- t.res$pval
        results <- rbind(results,  data.frame(sim=sim.name, iteration=i, stat=stat.name, n=n, pval=pval))
      }
    }
  }
})

saveRDS(results, file.path(opath, paste('discr_sims', '.rds', sep="")))
