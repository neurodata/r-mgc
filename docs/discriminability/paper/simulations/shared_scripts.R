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

## ------------------------------------------
# Simulations
## ------------------------------------------
sim_gmm <- function(mus, Sigmas, n, priors=NULL) {
  K <- dim(mus)[2]
  if (is.null(priors)) {
    priors <- rep(1/K, K)
  }
  ni <- rowSums(rmultinom(n, 1, prob=priors))
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

validator <- function(X, Y, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
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
