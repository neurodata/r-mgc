require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(energy)
require(tidyverse)
require(mltools)
require(data.table)
require(R.utils)
source('./data_xfms.R')

genomics.dat <- readRDS('../data/real/genomics_data.rds')
data <- list(CPM=genomics.dat$CPM, counts=genomics.dat$counts)
Subject <- genomics.dat$covariates$donor
genom.covs <- genomics.dat$covariates[, grepl("characteristic", names(genomics.dat$covariates))]
Covariates <- list(Sex=rowSums(genom.covs == "Sex: Male"),
                   Age=sapply(1:nrow(genom.covs), function(i) {
                     j <- as.numeric(which(sapply(genom.covs[i,], function(x) {
                       grepl('age: ', x) & !(grepl('passage', x))
                     })))
                     as.numeric(str_replace(genom.covs[i,j], 'age: ', ''))
                   }),
                   BMI=sapply(1:nrow(genom.covs), function(i) {
                     j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('bmi: ', as.character(x)))))
                     as.numeric(str_replace(as.character(genom.covs[i,j]), 'bmi: ', ''))
                   }),
                   Insulin.sens=sapply(1:nrow(genom.covs), function(i) {
                     j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('state: ', as.character(x)))))
                     as.numeric(as.character(str_replace(
                       as.character(genom.covs[i,j]), 'state: ', '')) == "Insulin sensitive")
                   }),
                   Race=as.matrix(one_hot(data.table(Race=factor(sapply(1:nrow(genom.covs), function(i) {
                     j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('race: ', as.character(x)))))
                     as.character(str_replace(as.character(genom.covs[i,j]), 'race: ', ''))
                   }))))))

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

stats <- list(Discr=discr.os, PICC=icc.os, I2C2=i2c2.os)

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
    withTimeout(do.call(experiment$alg, list(experiment$X, experiment$ID)), timeout=2000)
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
    return(list(X=X.xfm, dat.name=dat, task.name=task, Covariate=, xfm.name=xfm))
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

