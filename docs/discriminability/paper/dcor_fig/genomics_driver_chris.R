# the following can be run in the flashx docker container
# docker pull flashx/flashx
# docker run -ti --entrypoint /bin/bash flashx/flashx

# require(devtools)
# install_github('neurodata/r-mgc')
# install_github('neurodata/lol')
# install_github('muschellij2/I2C2')
# install.packages(c('energy', 'mltools', 'data.table',
#                    'R.utils'))

require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(energy)
require(mltools)
require(data.table)
require(R.utils)
require(FlashR)
require(combinat)

source('./data_xfms.R')
source('../simulations/shared_scripts.R')
dat.base <- '/genomics/files'

Sessions <- c("A", "B", "C", "D", "E", "F", "G", "H")
# manner in which a session will be aggregated by "summing" across wells
Agg.Sessions <- list(Small=list(Session1=c("A"), Session2=c("E")),
                     Medium=list(Session1=c("A", "B"), Session2=c("E", "F")),
                     Large=list(Session1=c("A", "B", "C", "D"), Session2=c("E", "F", "G", "H")),
                     Single=list(Session1=c("A"), Session2=c("B"), Session3=c("C"), Session4=c("D"),
                                 Session5=c("E"), Session6=c("F"), Session7=c("G"), Session8=c("H")))

Resolutions <- c("500kb", "50kb", "5MB", "chr", "amplicon")
# parse the data into usable format
if (!file.exists(file.path('/genomics', 'chris_parsed.rds'))) {
  parsed.dat <- lapply(Resolutions, function(resolution) {
    if (resolution == "amplicon") {
      ext <- ".snp"
    } else if (resolution == "chr") {
      ext <- ".transposed_fraction"
    } else {
      ext <- ".reads"
    }
    condition.dat <- lapply(c("cancer", "normal"), function(status) {
      file.names <- list.files(file.path(dat.base, status, resolution))
      # get non-aggregate files
      if (resolution == "amplicon") {
        file.names <- file.names[!grepl(paste0("faster1.snp"), file.names)]
      } else if (resolution == "chr") {
        file.names <- file.names[!grepl(paste0("faster1.transposed_fraction"), file.names)]
      } else {
        file.names <- file.names[!grepl(paste0("faster1.", resolution), file.names)]
      }
      # get individual names
      individuals <- unique(gsub('(_).*', '', file.names))
      # get only the individuals with all 8 sessions
      full_ses <- sapply(individuals, function(id) sum(grepl(id, file.names))) == 8
      individuals <- individuals[full_ses]
      agg.dat <- lapply(individuals, function(id) {
        # read the participants' 8 data points, and sum across A-D and E-H
        id.dat <- do.call(cbind, lapply(Sessions, function(ses) {
          if (resolution == "chr") {
            fname <- file.path(dat.base, status, resolution,
                               paste0(id, "_faster1.", ses, ext))
          } else if (resolution == "amplicon") {
            fname <- file.path(dat.base, status, resolution,
                               paste0(id, "_faster1.", ses, ext))
            return(as.numeric(read.csv(fname, sep="", header=FALSE)$V3))
          } else {
            fname <- file.path(dat.base, status, resolution,
                               paste0(id, "_faster1.", ses, ".", resolution, ext))
          }
          return(as.numeric(read.csv(fname, header=TRUE)[[1]]))
        }))
        colnames(id.dat) <- Sessions
        dat.agg=do.call(c, lapply(names(Agg.Sessions), function(aggregation) {
          Agg.Session=Agg.Sessions[[aggregation]]
          dat.thisagg <- do.call(cbind, lapply(Agg.Session, function(Session) {
            if (length(Session) == 1) {
              return(id.dat[,Session,drop=FALSE])
            } else if (resolution == "chr") {
              return(apply(id.dat[,Session,drop=FALSE], 1, mean))
            } else {
              return(apply(id.dat[,Session,drop=FALSE], 1, sum))
            }
          }))
          colnames(dat.thisagg) <- rep(id, ncol(dat.thisagg))

          return(list(Data=t(dat.thisagg), Sessions=names(Agg.Session),
                      Individual=rep(id, ncol(dat.thisagg)), Aggregation=rep(aggregation, ncol(dat.thisagg))))
        }))

        genomics.dat <- do.call(rbind, lapply(agg.dat, function(dat) dat$Data))
        session.dat <- do.call(c, lapply(agg.dat, function(dat) dat$Sessions))
        individuals.dat <- do.call(c, lapply(agg.dat, function(dat) dat$Individual))
        aggregation=do.call(c, lapply(agg.dat, function(dat) dat$Aggregation))
        return(list(Data=genomics.dat, Sessions=session.dat,
                    Individuals=individuals.dat, Cancer=rep(status, length(individuals.dat)),
                    Aggregation=aggregation))
        })
    })
    genomics.dat <- do.call(rbind, lapply(condition.dat, function(dat) dat$Data))
    session.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Sessions))
    individual.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Individual))
    status.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Cancer))
    aggregation.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Aggregation))
    # status is a 1 if have cancer, a 0 otherwise
    return(list(Data=genomics.dat, Sessions=session.dat, Individuals=individual.dat,
                Status=as.numeric(status.dat == "cancer"), Resolution = resolution,
                Aggregation=aggregation.dat))
  })
} else {
  parsed.dat <- readRDS(file.path('/genomics', 'chris_parsed.rds'))
}


flashx.embed <- function(X, A) {
  return(as.matrix(X %*% A))
}

flashx.decomp <- function(X, ncomp=0) {
  svdX <- fm.svd(X, nu=0, nv=ncomp)
  decomp=list(comp=svdX$v, val=svdX$d)
  return(decomp)
}

flashx.pca <- function(X, r, ...) {
  X <- fm.as.matrix(X)
  # mean center by the column mean
  d <- dim(X)[2]
  if (r > d) {
    stop(sprintf("The number of embedding dimensions, r=%d, must be lower than the number of native dimensions, d=%d", r, d))
  }
  # center the data
  Xc  <- sweep(X, 2, colMeans(X), '-')
  X.decomp <- flashx.decomp(Xc, ncomp=r)

  return(list(Xr = flashx.embed(X, X.decomp$comp), d=X.decomp$val, A=X.decomp$comp))
}


stats <- list(Discr=discr.os, PICC=icc.os, I2C2=i2c2.os, Kernel=ksamp.os, FPI=fpi.os,
              DISCO=disco.os, HSIC=hsic.os)#, manova.os)
xfms <- list(Raw=nofn.xfm, Rank=ptr.xfm, Log=log.xfm, Unit=unit.xfm, Center=center.xfm,
             UnitVar=unitvar.xfm, ZScore=zscore.xfm)

if (!file.exists('/genomics/genomics_prep.rds')) {
  experiments.base <- do.call(c, lapply(parsed.dat, function(dat.res) {
    print(dat.res$Resolution)
    result <- lapply(names(xfms), function(xfm) {
      print(xfm)
      X.xfm <- do.call(xfms[[xfm]], list(dat.res$Data))
      print(dim(X.xfm))
      X.fm <- fm.as.matrix(X.xfm)
      DX <- as.matrix(fm.inner.prod(X.fm, t(X.fm), fm.bo.euclidean, fm.bo.add))
      RX <- as.matrix(cov(t(X.fm)))
      Xr <- as.matrix(flashx.pca(X.fm, 1)$Xr)
      rm(X.fm)
      gc()
      return(list(X=X.xfm, DX=DX, Xr=Xr, RX=RX, Individuals=dat.res$Individuals, Cancer=dat.res$Status,
                  xfm.name=xfm, Resolution=dat.res$Resolution, Sessions=dat.res$Sessions))
    })
    gc()
    return(result)
  }))

  saveRDS(experiments.base, '/genomics/genomics_prep.rds')
} else {
  experiments.base <- readRDS('/genomics/genomics_prep.rds')
}


results.reference <- do.call(rbind, mclapply(experiments.base, function(experiment) {
  print(sprintf("Resolution=%s, XFM=%s", experiment$Resolution,
                experiment$xfm.name))
  do.call(rbind, lapply(names(stats), function(stat.name) {
    tryCatch({
      if (stat.name %in% c("Discr", "Kernel", "HSIC")) {
        X.dat = experiment$DX
      } else if (stat.name == "PICC") {
        X.dat = experiment$Xr
      } else if (stat.name == "FPI") {
        X.dat = experiment$RX
      } else {
        X.dat = experiment$X
      }
      if (stat.name == "FPI") {
        stat.res=do.call(stats[[stat.name]], list(X=X.dat, Y=experiment$Individual, Z=experiment$Sessions, is.sim_or_dist=TRUE, is.sim=TRUE))
      } else {
        stat.res=do.call(stats[[stat.name]], list(X.dat, experiment$Individual))
      }
    }, error=function(e) {
      print(sprintf("Resolution=%s, XFM=%s, Stat=%s, ERROR=%s", experiment$Resolution,
                    experiment$xfm.name, stat.name, e))
      return(NULL)
      })
    return(data.frame(Resolution=experiment$Resolution, xfm=experiment$xfm.name,
                      Algorithm=stat.name, Statistic=stat.res))
  }))
}, mc.cores=detectCores()))

results.effect <- do.call(rbind, mclapply(experiments.base, function(experiment) {
  print(sprintf("Resolution=%s, XFM=%s", experiment$Resolution, experiment$xfm.name))
  stat=tryCatch({
    ksamp.test(experiment$DX, experiment$Cancer, nrep=1000L)
  }, error=function(e) {
    print(sprintf("Data=%s, XFM=%s, ERROR=%s", experiment$dat.name, experiment$xfm.name, e))
    return(NULL)
  })
  if (!is.null(stat)) {
    return(data.frame(Resolution=experiment$Resolution, xfm=experiment$xfm.name, #Aggregation=experiment$Aggregation,
                      Algorithm="DCorr", statistic=stat$statistic, pvalue=stat$pvalue))
  } else {
    return(NULL)
  }
}, mc.cores=detectCores()))

saveRDS(list(Reference=results.reference, Effect=results.effect), '../data/real/genomics_cancer.rds')
