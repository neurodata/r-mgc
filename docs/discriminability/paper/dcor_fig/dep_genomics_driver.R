require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(igraph)
require(fmriutils)
require(reshape2)
require(stringr)
require(FNN)
require(Metrics)
require(randomForest)
require(rARPACK)
require(energy)



stats <- list(discr.os, anova.os, icc.os, i2c2.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2")

graph.xfms <- list(nofn, ptr, log.xfm)
names(graph.xfms) <- c("N", "P", "L")

genomics.dat <- readRDS('../data/real/genomics_data.rds')
data <- list(CPM=genomics.dat$CPM, counts=genomics.dat$counts)
Subject <- genomics.dat$covariates$donor
Y <- genomics.dat$covariates$

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
    return(log2(x + min(x)/2))
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

stats <- list(Discr=discr.os, ANOVA=anova.os, PICC=icc.os, I2C2=i2c2.os)

xfms <- list(Raw=nofn.xfm, Rank=ptr.xfm, Log=log.xfm, Unit=unit.xfm, Center=center.xfm,
             UnitVar=unitvar.xfm, ZScore=zscore.xfm)

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  do.call(rbind, lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], X)
    lapply(names(stats), function(stat) {
      return(X=X.xfm, dat.name=dat, ID=Subject,
             xfm.name=xfm, alg=algs[[stat]], alg.name=stat)
    })
  }))
}))

mclapply(experiments, function(experiment) {
  X <- do.call(experiment$xfm, experiment$X)
  stat <- do.call(experiment$alg, list(X, experiment$ID))
  return(data.frame(Data=experiment$dat.name, xfm=experiment$xfm.name, Algorithm=experiment$alg.name,
                    Statistic=stat))
}, mc.cores=detectCores() - 1)

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  lapply(names(xfms), function(xfm) {
  X.xfm <- do.call(xfms[[xfm]], X)
    return(X=X.xfm, dat.name=dat, Sex=Sex, xfm.name=xfm)
  })
}))

mclapply(experiments, function(experiment) {
  do.call(dcor, list(experiment$X, experiment$Sex))
}, mc.cores=detectCores() - 1)
