##----------------------------------------------
# One Sample Tests
##----------------------------------------------
require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
no_cores = detectCores() - 1


# one-way ICC
icc.os <- function(X, y) {
  data <- data.frame(x=lol.project.pca(X, 1), y=y)
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
test.one_sample <- function(X, Y, is.dist=FALSE, dist.xfm=discr.distance, dist.params=list(method='euclidean'),
                            dist.return=NULL, remove.isolates=TRUE, nperm=100, no_cores=1) {
  
  validated <- discr.validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
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
               pval=(sum(nr[[stat.name]]>tr[[stat.name]]) + 1)/(nperm + 1))
  })))
}
