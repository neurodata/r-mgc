
nofn.xfm <- function(x, ...) {
  return(x)
}

ptr.xfm <- function(X, ...) {
  ptr.col <- function(x) {
    nz <- x[x != 0]
    r <- rank(nz)*2/(length(nz) + 1)
    x[x != 0] <- r
    x <- (x - min(x))/(max(x) - min(x))
    return(x)
  }
  return(apply(X, 2, ptr.col))
}

log.xfm <- function(X, ...) {
  log.col <- function(x) {
    return(log2(x + min(x[x != 0])/2))
  }
  return(apply(X, 2, log.col))
}

unit.xfm <- function(X, ...) {
  unit.col <- function(x) {
    return((x - min(x))/(max(x)-min(x)))
  }
  return(apply(X, 2, unit.col))
}

center.xfm <- function(X, ...) {
  center.col <- function(x) {
    return(x - mean(x))
  }
  return(apply(X, 2, center.col))
}

unitvar.xfm <- function(X, ...) {
  unitvar.col <- function(x) {
    return(x/sd(x))
  }
  return(apply(X, 2, unitvar.col))
}

zscore.xfm <- function(X, ...) {
  zsc.col <- function(x) {
    x.m <- x - mean(x)
    return((x.m)/sd(x.m))
  }
  return(apply(X, 2, zsc.col))
}