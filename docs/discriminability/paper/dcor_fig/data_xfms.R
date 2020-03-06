
nofn.xfm <- function(x, ...) {
  return(x)
}

ptr.xfm <- function(X, ...) {
  ptr.row <- function(x) {
    nz <- x[x != 0]
    r <- rank(nz)*2/(length(nz) + 1)
    x[x != 0] <- r
    x <- (x - min(x))/(max(x) - min(x))
    return(x)
  }
  return(t(apply(X, 1, ptr.row)))
}

log.xfm <- function(X, ...) {
  log.row <- function(x) {
    if (length(unique(x)) != 1) {
      return(log2(x + min(x[x != 0])/2))
    } else {
      return(rep(0, length(x)))
    }
  }
  return(t(apply(X, 1, log.row)))
}

unit.xfm <- function(X, ...) {
  unit.row <- function(x) {
    if (length(unique(x)) != 1) {
      return((x - min(x))/(max(x)-min(x)))
    } else {
      return(rep(0, length(x)))
    }
  }
  return(t(apply(X, 1, unit.row)))
}

center.xfm <- function(X, ...) {
  center.row <- function(x) {
    return(x - mean(x))
  }
  return(t(apply(X, 1, center.row)))
}

unitvar.xfm <- function(X, ...) {
  unitvar.row <- function(x) {
    if (length(unique(x)) != 1) {
      return(x/sd(x))
    } else {
      return(rep(0, length(x)))
    }
  }
  return(t(apply(X, 1, unitvar.row)))
}

zscore.xfm <- function(X, ...) {
  zsc.row <- function(x) {
    if (length(unique(x)) != 1) {
      x.m <- x - mean(x)
      return((x.m)/sd(x.m))
    } else {
      return(rep(0, length(x)))
    }
  }
  return(t(apply(X, 1, zsc.row)))
}
