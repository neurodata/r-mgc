require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(energy)
require(tidyverse)


genomics.dat <- readRDS('../data/real/genomics_data.rds')
data <- list(CPM=genomics.dat$CPM, counts=genomics.dat$counts)
Subject <- genomics.dat$covariates$donor
genom.covs <- genomics.dat$covariates[, grepl("characteristic", names(genomics.dat$covariates))]
Sex <- rowSums(genom.covs == "Sex: Male")

# one-way ICC
icc.os <- function(X, y) {
  x <- lol.project.pca(X, r=1)$Xr
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

# one-sample MANOVA
manova.os <- function(X, Y) {
  fit <- manova(X ~ Y)
  return(summary(fit)$stats["Y", "approx F"])
}

# I2C2 wrapper
i2c2.os <- function(X, Y) {
  return(I2C2.original(y=X, id=Y, visit=rep(1, length(Y)), twoway=FALSE)$lambda)
}

discr.os <- function(X, Y) {
  return(discr.stat(X, Y)$discr)
}

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

stats <- list(Discr=discr.os)

xfms <- list(Raw=nofn.xfm, Rank=ptr.xfm, Log=log.xfm, Unit=unit.xfm, Center=center.xfm,
             UnitVar=unitvar.xfm, ZScore=zscore.xfm)

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  do.call(rbind, lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], list(X))
    lapply(names(stats), function(stat) {
      return(list(X=X.xfm, dat.name=dat, ID=Subject,
             xfm.name=xfm, alg=stats[[stat]], alg.name=stat))
    })
  }))
}))

result.stat <- do.call(rbind, lapply(experiments, function(experiment) {
  print(sprintf("Data=%s, XFM=%s, Alg=%s", experiment$dat.name, experiment$xfm.name, experiment$alg.name))
  stat <- tryCatch({
    withTimeout(do.call(experiment$alg, list(experiment$X, experiment$ID)), timeout=1200)
  }, error=function(e) {
    print(e)
    return(NULL)
  })
  if (!is.null(stat)) {
    return(data.frame(Data=experiment$dat.name, xfm=experiment$xfm.name, Algorithm=experiment$alg.name,
                      Statistic=stat))
  } else {
    return(NULL)
  }
}))

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], list(X))
    return(list(X=X.xfm, dat.name=dat, Sex=Sex, xfm.name=xfm))
  })
}))

result.dcor <- do.call(rbind, mclapply(experiments, function(experiment) {
  stat=tryCatch({
    do.call(dcor, list(experiment$X, experiment$Sex))
  }, error=function(e) {
    print(sprintf("Data=%s, XFM=%s", experiment$dat.name, experiment$xfm.name))
    return(NULL)
  })
  if (!is.null(stat)) {
    return(data.frame(Data=experiment$dat.name, xfm=experiment$xfm.name, Algorithm="DCor",
                      Statistic=stat))
  } else {
    return(NULL)
  }
}, mc.cores=detectCores() - 1))

