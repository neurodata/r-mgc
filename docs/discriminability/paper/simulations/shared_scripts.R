require(mgc)
require(lolR)
require(MASS)
require(mvtnorm)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(lolR)
require(abind)
require(emdbook)
require(gtools)
require(kernelPSI)

# one-way ICC
icc.os <- function(X, Y, ...) {
  if (length(dim(X)) == 2) {
    if (dim(X)[2] > 1) {
      X=lol.project.pca(X, 1)$Xr
    }
  }
  data <- data.frame(X=X, Y=Y, ...)
  fit <- anova(aov(X ~ Y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  a <- length(unique(Y))
  tmp.outj <- as.numeric(aggregate(X ~ Y, data=data, FUN = length)$X)
  k <- (1/(a - 1)) * (sum(tmp.outj) - (sum(tmp.outj^2)/sum(tmp.outj)))
  var.a <- (MSa - MSw)/k
  r <- var.a/(var.w + var.a)
  return(r)
}

# I2C2 wrapper
i2c2.os <- function(X, Y, ...) {
  return(I2C2.original(y=X, id=Y, visit=rep(1, length(Y)), twoway=FALSE)$lambda)
}


# average fingerprint index
# X is a nxd data matrix for n samples
# Y is a n vector of individual labels
# Z is a n vector of scan sessions
fpi.os <- function(X, Y, Z, ...) {
  DX <- cor(t(X))  # compute pairwise correlation
  # compute the ranks per scan
  RX <- apply(DX, 1, function(x) rank(-x, ties.method="first"))
  individuals <- unique(Y)
  sessions <- unique(Z)
  fpis <- unlist(sapply(1:length(individuals), function(i) {
    individual <- individuals[i]
    # find the indices associated with current individual
    idx.i <- which(Y == individual)
    # get the session indices associated with current individual
    sessions.i <- Z[idx.i]
    # compute all permutations of the sessions available for individual
    ses.perm <- permutations(length(sessions.i), 2)
    # for each (session i, session j) pair, compute fingerprint index
    # for this individual
    sapply(1:dim(ses.perm)[1], function(q) {
      # get the database sessions
      ses.idx <- which(Z == sessions.i[ses.perm[q,2]])
      # if there are no other items with this particular session, skip
      if (length(ses.idx) == 0) {
        return(NaN)
      }
      # restrict RX to the database session, and compute the index
      # of the minimum
      min.idx <- which.min(RX[idx.i[ses.perm[q,1]], ses.idx])
      # compare the subject associated with the minimum to the
      # current subject; if so, return 1, else 0
      ifelse(individual %in% Y[ses.idx][min.idx], return(1), return(0))
    })
  }))
  return(mean(fpis, na.rm=TRUE))
}

discr.os <- function(X, Y, is.dist=TRUE, ...) {
  return(discr.stat(X, Y, is.dist=is.dist)$discr)
}

disco.os <- function(X, Y, is.dist=TRUE, ...) {
  if (is.dist) {
    DX <- X
  } else {
    DX <- mgc.distance(X, method="euclidean")
  }
  as.numeric(disco(X, factor(Y), R=0, method="disco")$statistic)
}

mmd.os <- function(X, Y, is.dist=FALSE) {
  if (is.dist) {
    DX <- X
  } else {
    DX <- mgc.distance(X, method="euclidean")
  }
  DY <- mgc.distance(Y, method="ohe")
  return(HSIC(DX, DY))
}

## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n, priors=NULL) {
  K <- dim(mus)[2]
  if (is.null(priors)) {
    ni <- rep(floor(n/K), K)
  } else {
    ni <- rowSums(rmultinom(n, 1, prob=priors))
  }
  X <- do.call(rbind, lapply(1:K, function(k) matrix(mvrnorm(n=ni[k], mus[,k], Sigmas[,,k]), nrow=ni[k])))
  Y <- do.call(c, lapply(1:K, function(k) rep(k, ni[k])))
  return(list(X=X, Y=Y))
}

sim_gmm_match <- function(mus, Sigmas, Ys) {
  K <- length(sort(unique(Ys)))
  ni <- sapply(sort(unique(Ys)), function(y) sum(Ys == y))
  X <- do.call(rbind, lapply(1:K, function(k) matrix(mvrnorm(n=ni[k], mus[,k], Sigmas[,,k]), nrow=ni[k])))
  return(list(X=X, Y=Ys))
}

validator <- function(X, Y, Z, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                      dist.return=NULL, remove.isolates=TRUE) {

  DX <- mgc:::mgc.dist.validator(X, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return)

  # validate Y
  if (!is.vector(Y)) {
    tryCatch({
      Y <- as.vector(Y)
    }, error=function(e) stop("You have not passed an object Y that can be coerced to a [n] vector."))
  }

  if (nrow(DX) != length(Y)) {
    stop("Your distance matrix and your ids vector do not have the same number of samples, after applying distance function.")
  }

  # remove isolated subjects if requested.
  if (remove.isolates) {
    purged <- mgc:::remove.isolates(DX, Y, is.dist=TRUE)
    DX <- purged$X; Y <- purged$Y
    purged.X <- mgc:::remove.isolates(X, Y, is.dist=FALSE)
    X <- purged$X; Y.p <- purged$Y
    if (!all(Y == Y.p)) {
      stop("Something went amuck.")
    }
  }


  if (length(unique(Y)) <= 1) {
    stop("You have passed a vector containing only a single unique sample id.")
  }

  return(list(D=DX, X=X, Y=Y))
}
