
# Parallelize Stuff
#==================
require(parallel)
require(mgc)
require(plyr)

no_cores = detectCores() - 1

# Execute Simulations
#====================

# Setup Algorithms
#=================
sim_types <- list(discr.sims.linear, discr.sims.log, discr.sims.spread, discr.sims.radial)
sim_names <- c("Linear", "Logarithmic", "Spread", "Radial")
names(sim_types) <- sim_names

Ks <- c(2, 10, 10)
class.equals <- c(TRUE, TRUE, FALSE)
nrep <- 20
nperm <- 100
ds <- c(1, 100)
nmin <- 30
nmax <- 1000
length.out=15

log.seq <- function(a=20, b=1000, base=2, length.out=20) {
  return(round(base^seq(log(a, base=base), log(b, base=base), length.out=length.out)))
}
algs <- list(discr.stat, i2c2)
names(algs) <- c("Discriminability", "I2C2")

fold_rep <- data.frame(n=c(), p=c(), simname=c(), iteration=c(), K=c(), class.equal=c(), alg_name=c())
counter <- 1

for (i in 1:length(sim_types)) {
  n <- rep(log.seq(a=nmin, b=nmax, base=2, length.out=length.out), nrep*length(ds))
  iteration <- rep(rep(seq(1, length.out), each=nrep), length(ds))
  d <- rep(ds, each=nrep*length.out)
  for (k in 1:length(Ks)) {
    fold_rep <- rbind(fold_rep, data.frame(n=n, p=d, simname=sim_names[[i]], iteration=iteration, K=Ks[k],
                                           class.equal=class.equals[k]))
  }
}

results <- mclapply(fold_rep, function(sim) {
  tryCatch({
    simdat <- do.call(sim_types[[sim$simname]], list(n=sim$n, d=sim$p, K=sim$K, class.equal=sim$class.equal))
    perm.test <- discr.test.one_sample(simdat$X, simdat$Y, nperm=nperm)
    res <- data.frame(n=sim$n, p=sim$p, simname=sim$simname, iteration=sim$iteration, K=sim$K,
                      class.equal=sim$class.equal, tstat=perm.test$srel, pval=perm.test$pval)
    return(res)
  }, error=function(e) lhat <- NaN)
}, mc.cores=no_cores)

resultso <- do.call(rbind, results)
saveRDS(resultso, file.path(opath, paste('discr_sims.rds')))
