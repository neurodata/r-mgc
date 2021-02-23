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
source('../simulations/shared_scripts.R')

genomics.dat <- readRDS('../data/real/genomics_data.rds')
data <- list(CPM=genomics.dat$CPM, counts=genomics.dat$counts)
Subject <- genomics.dat$covariates$donor
Session <- as.numeric(genomics.dat$covariates$lineRep)
genom.covs <- genomics.dat$covariates[, grepl("characteristic", names(genomics.dat$covariates))]
# Covariates <- list(Sex=rowSums(genom.covs == "Sex: Male"),
#                    Age=sapply(1:nrow(genom.covs), function(i) {
#                      j <- as.numeric(which(sapply(genom.covs[i,], function(x) {
#                        grepl('age: ', x) & !(grepl('passage', x))
#                      })))
#                      as.numeric(str_replace(genom.covs[i,j], 'age: ', ''))
#                    }),
#                    BMI=sapply(1:nrow(genom.covs), function(i) {
#                      j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('bmi: ', as.character(x)))))
#                      as.numeric(str_replace(as.character(genom.covs[i,j]), 'bmi: ', ''))
#                    }),
#                    Insulin.sens=sapply(1:nrow(genom.covs), function(i) {
#                      j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('state: ', as.character(x)))))
#                      as.numeric(as.character(str_replace(
#                        as.character(genom.covs[i,j]), 'state: ', '')) == "Insulin sensitive")
#                    }),
#                    Race=as.matrix(one_hot(data.table(Race=factor(sapply(1:nrow(genom.covs), function(i) {
#                      j <- as.numeric(which(sapply(genom.covs[i,], function(x) grepl('race: ', as.character(x)))))
#                      as.character(str_replace(as.character(genom.covs[i,j]), 'race: ', ''))
#                    }))))))
Sex=rowSums(genom.covs == "Sex: Male")

stats <- list(Discr=discr.os, PICC=icc.os, I2C2=i2c2.os, Kernel=ksamp.os, FPI=fpi.os,
              DISCO=disco.os, HSIC=hsic.os)#, manova.os)

xfms <- list(Raw=nofn.xfm, Rank=ptr.xfm, Log=log.xfm, Unit=unit.xfm, Center=center.xfm,
             UnitVar=unitvar.xfm, ZScore=zscore.xfm)

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  do.call(rbind, lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], list(X))
    lapply(names(stats), function(stat) {
      return(list(X=X.xfm, dat.name=dat, ID=Subject, Session=Session,
                  xfm.name=xfm, alg=stats[[stat]], alg.name=stat))
    })
  }))
}))

result.stat <- do.call(rbind, mclapply(experiments, function(experiment) {
  print(sprintf("Data=%s, XFM=%s, Alg=%s", experiment$dat.name, experiment$xfm.name, experiment$alg.name))
  stat <- tryCatch({
    if (experiment$alg.name != "FPI") {
      withTimeout(do.call(experiment$alg, list(experiment$X, experiment$ID, is.dist=FALSE)), timeout=2000)
    } else {
      withTimeout(do.call(experiment$alg, list(experiment$X, experiment$ID, experiment$Session, is.dist=FALSE)), timeout=5000)
    }
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
  return(stat)
}, mc.cores=detectCores() - 1))

experiments <- do.call(rbind, lapply(names(data), function(dat) {
  X <- t(data[[dat]])
  lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], list(X))
    return(list(X=X.xfm, dat.name=dat, Sex=Sex, xfm.name=xfm))
  })
}))

result.dcor <- do.call(rbind, mclapply(experiments, function(experiment) {
  stat=tryCatch({
    do.call(ksamp.test, list(experiment$X, experiment$Sex))
  }, error=function(e) {
    print(sprintf("Data=%s, XFM=%s", experiment$dat.name, experiment$xfm.name))
    return(NULL)
  })
  if (!is.null(stat)) {
    return(data.frame(Data=experiment$dat.name, xfm=experiment$xfm.name, Algorithm="DCorr",
                      statistic=stat$statistic, pvalue=stat$pvalue))
  } else {
    return(NULL)
  }
  return(stat)
}, mc.cores=detectCores() - 1))

saveRDS(list(Reference=result.stat, Effect=result.dcor), '../data/real/genomics_sex.rds')
