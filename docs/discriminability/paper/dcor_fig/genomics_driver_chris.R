# the following can be run in the flashx docker container
# docker pull flashx/flashx
# docker run -ti --entrypoint /bin/bash flashx/flashx

require(devtools)
install_github('neurodata/r-mgc')
install_github('neurodata/lol')
install_github('muschellij2/I2C2')
install.packages(c('energy', 'mltools', 'data.table',
                   'R.utils'))


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

source('./data_xfms.R')
dat.base <- '/genomics/files'

Sessions <- c("A", "B", "C", "D", "E", "F", "G", "H")
# manner in which a session will be aggregated by "summing" across wells
Agg.Sessions <- list(Session1=c("A"), 
                     Session2=c("B"),
                     Session3=c("C"),
                     Session4=c("D"),
                     Session5=c("E"),
                     Session6=c("F"),
                     Session7=c("G"),
                     Session8=c("H"))
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
        if (!all(sapply(Agg.Sessions, function(agg) length(agg) == 1))) {
          id.dat <- do.call(cbind, lapply(Agg.Sessions, function(Session) {
            return(apply(id.dat[,Session,drop=FALSE], 1, sum))
          }))
        }
        colnames(id.dat) <- rep(id, ncol(id.dat))
        
        return(list(Data=t(id.dat), Sessions=names(Agg.Sessions), Individual=id))
      })
      genomics.dat <- do.call(rbind, lapply(agg.dat, function(dat) dat$Data))
      session.dat <- do.call(c, lapply(agg.dat, function(dat) dat$Sessions))
      individuals.dat <- do.call(c, lapply(agg.dat, function(dat) rep(dat$Individual, length(dat$Sessions))))
      return(list(Data=genomics.dat, Sessions=session.dat, 
                  Individuals=individuals.dat, Status=rep(status, length(individuals.dat))))
    })
    genomics.dat <- do.call(rbind, lapply(condition.dat, function(dat) dat$Data))
    session.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Sessions))
    individual.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Individual))
    status.dat <- do.call(c, lapply(condition.dat, function(dat) dat$Status))
    # status is a 1 if have cancer, a 0 otherwise
    return(list(Data=genomics.dat, Sessions=session.dat, Individuals=individual.dat, 
                Status=as.numeric(status.dat == "cancer"), Resolution = resolution))
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


# one-way ICC
icc.os <- function(X, y) {
  data <- data.frame(x=X, y=y)
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

discr.os <- function(X, Y, is.dist=FALSE) {
  return(discr.stat(X, Y, is.dist=is.dist)$discr)
}

disco.os <- function(X, Y, is.dist=FALSE) {
  if (is.dist) {
    DX <- X
  } else {
    DX <- mgc.distance(X, method="euclidean")
  }
  disco(X, factor(Y))
}

stats <- list(SimilRR=discr.os, PICC=icc.os, I2C2=i2c2.os, DISCO=disco.os)
xfms <- list(Raw=nofn.xfm, Rank=ptr.xfm, Log=log.xfm, Unit=unit.xfm, Center=center.xfm,
             UnitVar=unitvar.xfm, ZScore=zscore.xfm)

experiments.base <- do.call(c, lapply(parsed.dat, function(dat.res) {
  result <- lapply(names(xfms), function(xfm) {
    X.xfm <- do.call(xfms[[xfm]], list(dat.res$Data))
    X.fm <- fm.as.matrix(X.xfm)
    DX <- as.matrix(fm.inner.prod(X.fm, t(X.fm), fm.bo.euclidean, fm.bo.add))
    Xr <- as.matrix(flashx.pca(X.xfm, 1)$Xr)
    rm(X.fm)
    return(list(X=X.xfm, DX=DX, Xr=Xr, Individuals=dat.res$Individuals, Status=dat.res$Status,
                xfm.name=xfm, Resolution=dat.res$Resolution))
  })
  gc()
  return(result)
}))

saveRDS(experiments.base, '../data/real/genomics_prep.rds')