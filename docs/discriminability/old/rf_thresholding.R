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
# load/source MASE code
mase.path <- './mase/R/'
mase.files <- list.files(mase.path)
mase.files <- mase.files[mase.files %in% c("mase.R", "omnibus-embedding.R", "getElbows.R")]
sapply(mase.files, function(x) source(file.path(mase.path, x)))

#fmri.path <- '/mnt/nfs2/MR/cpac_3-9-2/'
#pheno.path <- '/mnt/nfs2/MR/all_mr/phenotypic/'
fmri.path <- '/cis/project/ndmg/eric/discriminability/cpac_3-9-2/'
pheno.path <- '/cis/project/ndmg/eric/discriminability/phenotypic/'
#fmri.path <- '/data/cpac_3-9-2/'
#pheno.path <- '/data/all_mr/phenotypic/'
opath <- './data/real/'
no_cores <- parallel::detectCores() - 2

# one-way anova
anova.os <- function(X, y) {
  x <- lol.project.pca(X, r=1)$Xr
  data <- data.frame(x=x, y=y)
  fit <- anova(aov(x ~ y, data=data))
  MSa <- fit$"Mean Sq"[1]
  MSw <- var.w <- fit$"Mean Sq"[2]
  f = fit[["F value"]][1]
  p = fit[["Pr(>F)"]][1]
  return(f)
}

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
manova.onesample.driver <- function(X, Y) {
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

cpac.open_graphs <- function(fnames, dataset_id="", atlas_id="",
                             fmt='elist', verbose=FALSE, rtype='list', flatten=FALSE,
                             rem.diag=TRUE, sub_pos=2, ses_pos=4) {
  if (! (fmt %in% c('adj', 'elist', 'graphml'))) {
    stop('You have passed an invalid format type. Options are: [\'adj\', \'elist\', and \'graphml\'].')
  }

  if (fmt == 'elist') {
    fmt = 'ncol'; ext = "ssv"
  } else if (fmt == "graphml") {
    fmt = "graphml"; ext = "graphml"
  } else if (fmt == "adj") {
    fmt = "adj"; ext="adj"
  }

  if (is.character(fnames)) {
    fnames <- list.files(fnames, pattern=paste('\\.', ext, sep=""), full.names=TRUE)
  }

  if (! (rtype %in% c('list', 'array'))) {
    stop('You have passed an invalid return type. Options are: [\'list\', \'array\'].')
  }

  print(sprintf("opening graphs for %s dataset and %s parcellation atlas...", dataset_id, atlas_id))
  subjects <- vector("character", length(fnames))
  sessions <- vector("character", length(fnames))
  tasks <- vector("character", length(fnames))
  gr <- list()

  vertices <- c()

  # so that we don't get any annoying errors if particular vertices are empty
  if (fmt != "adj") {
    for (i in 1:length(fnames)) {
      tgr <- igraph::read_graph(fnames[i], format=fmt) # read the graph from the filename
      vertices <- base::union(vertices, as.numeric(V(tgr)$name))
    }
  }

  vertices <- sort(vertices)
  counter <- 1
  for (i in 1:length(fnames)) {
    basename <- basename(fnames[i])     # the base name of the file
    if (verbose) {
      print(paste('Loading', basename, '...'))
    }
    tgr <- tryCatch({
      igraph::read_graph(fnames[i], format=fmt, predef=vertices) # read the graph from the filename, ordering by the vertices we found previously
    }, error = function(e) {
      return(NaN)
    })

    if (is.igraph(tgr)) {
      tgr <- get.adjacency(tgr, type="both", attr="weight", sparse=FALSE) # convert graph to adjacency matrix
      tgr[is.nan(tgr)] <- 0  # missing entries substituted with 0s
      if (rem.diag) {
        diag(tgr) <- 0
      }
      gr[[basename]] <-t(tgr)
      str_names <- strsplit(basename, '_')[[1]]
      subjects[counter] <- str_names[sub_pos]
      sessions[counter] <- str_names[ses_pos]
      counter <- counter + 1
    }
  }

  dataset <- rep(dataset_id, counter - 1)
  atlas <- rep(atlas_id, counter - 1)
  subjects <- subjects[1:counter - 1]
  sessions <- sessions[1:counter - 1]

  if (rtype == 'array') {
    aro <- fmriu.list2array(gr, flatten=flatten)
    gr <- aro$array
    dataset <- dataset[aro$incl_ar]
    atlas <- atlas[aro$incl_ar]
    subjects <- subjects[aro$incl_ar]
    sessions <- sessions[aro$incl_ar]
  }
  return(list(graphs=gr, dataset=dataset, atlas=atlas, subjects=subjects,
              sessions=sessions))
}

fmriu.list2array <- function(list_in, flatten=FALSE) {
  nroi <- max(sapply(list_in, function(graph) dim(graph)[1]))
  nsub <- length(list_in)
  array_out <- array(NaN, dim=c(nsub, nroi, nroi))
  subnames <- names(list_in)
  incl_ar <- logical(nsub)
  for (i in 1:nsub) {
    if (isTRUE(all.equal(dim(list_in[[i]]), c(nroi, nroi)))) {
      array_out[i,,] <-list_in[[i]]
      incl_ar[i] <- TRUE
    }
  }
  array_out <- array_out[incl_ar,,]
  subnames <- subnames[incl_ar]
  if (flatten) {
    dimar <- dim(array_out)
    dim(array_out) <- c(dimar[1], dimar[2]*dimar[3])
  }
  return(list(array=array_out, incl_ar=incl_ar, names=subnames))
}

dsets <- list.dirs(path=fmri.path, recursive=FALSE)

# atlas_opts <- c("C", "D")
# names(atlas_opts) <- c("cc2", "des")
atlas_opts <- c("D")
names(atlas_opts) <- c("des")


trng <- seq(0, 1, by=0.025)

# run all datasets at all parcel resolutions
experiments <- do.call(c, lapply(dsets, function(dset) {
  dset_name = basename(dset)
  if (dset_name %in% c("MPG1", "BNU3")) {
    return(NULL)
  }
  if (grepl("NKI24", dset_name)) {
    if (!grepl("std2500", dset_name)) {
      return(NULL)
    }
    dset.key <- "NKI1"
    sub.pos <- 3
  } else if (grepl("KKI", dset_name)) {
    dset.key <- dset_name
    sub.pos <- 1
  } else {
    dset.key <- dset_name
    sub.pos <- 2
  }
  do.call(c, res.thresh <- lapply(trng, function(thr) {
    lapply(names(atlas_opts), function(atlas) {
      list(Dataset=dset.key, Parcellation=atlas_opts[atlas],
           dat.path=file.path(fmri.path, dset_name, "graphs", sprintf("FSL_nff_nsc_gsr_%s", atlas)),
           pheno.path=file.path(pheno.path, paste(dset.key, "_phenotypic_data.csv", sep="")),
           sub.pos = sub.pos, thr=thr)
    })
  }))
}))

stats <- list(discr.os, anova.os, icc.os, i2c2.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2")

rf.res.path <- file.path(opath, 'rf_results')
dir.create(opath)
dir.create(rf.res.path)

#================
#
# fMRI Driver
#
#================

# range of thresholds to try
# multicore apply over dataset
rf.results <- mclapply(experiments, function(exp) {
  graphs <- cpac.open_graphs(exp$dat.path, dataset_id=exp$Dataset,
                             atlas_id=exp$Parcellation, sub_pos = exp$sub.pos, flatten=FALSE)

  print(sprintf("Dataset: %s, thr=%.3f", exp$Dataset, exp$thr))

  # threshold the graphs
  graphs.bin <- lapply(graphs$graphs, function(graph) {
    tmp <- graph
    tmp[tmp < exp$thr] <- 0
    tmp[tmp >= exp$thr] <- 1
    return(tmp)
  })
  # flatten for statistics
  flat.gr <- fmriu.list2array(graphs.bin, flatten=TRUE)

  # run statistics
  res <- do.call(rbind, lapply(names(stats), function(stat) {
    tryCatch({
      return(data.frame(Dataset=exp$Dataset, thresh=exp$thr, alg=stat,
                        nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                        nroi=sqrt(dim(flat.gr$array)[2]), nsub=length(unique(graphs$subjects)),
                        stat=do.call(stats[[stat]], list(flat.gr$array, graphs$subjects))))
    }, error=function(e) {return(NULL)})
  }))

  graphs.embedded <- list(
    raw=t(simplify2array(lapply(graphs.bin, function(x) as.vector(x)))),
    mase=t(simplify2array(lapply(mase(graphs.bin)$R, function(x) as.vector(x))))
  )
  graphs.embedded$dist <- g.ase(as.matrix(dist(graphs.embedded$mase)))$X

  pheno.dat <- read.csv(exp$pheno.path)
  pheno.dat$AGE_AT_SCAN_1 <- as.numeric(as.character(pheno.dat$AGE_AT_SCAN_1))
  pheno.dat <- pheno.dat[!duplicated(pheno.dat$SUBID),]
  pheno.dat <- pheno.dat[, c("SUBID", "AGE_AT_SCAN_1", "SEX")]
  pheno.scans <- pheno.dat[sapply(as.numeric(graphs$subjects), function(x) which(x == pheno.dat$SUBID)),]
  if (exp$Dataset == "KKI2009") {
    pheno.scans$SEX <- as.factor((pheno.scans$SEX == "M") + 1)
  }

  task.res <- do.call(rbind, lapply(names(graphs.embedded), function(embed) {
    embed.graphs <- graphs.embedded[[embed]]
    tryCatch({

      # aggregate results across all subjects; report RMSE at current r
      age.res <- do.call(rbind, lapply(unique(graphs$subjects), function(sub) {
        training.set <- which(graphs$subjects != sub)  # hold out same-subjects from training set
        testing.set <- which(graphs$subjects == sub)  # validate over all scans for this subject
        # predict for held-out subject
        trained.age.rf <- randomForest(embed.graphs[training.set,], y=as.numeric(pheno.scans$AGE_AT_SCAN_1[training.set]))
        preds.age.rf <- predict(trained.age.rf, embed.graphs[testing.set,])
        return(data.frame(true=pheno.scans$AGE_AT_SCAN_1[testing.set],
                          pred=preds.age.rf, subject=sub))
      }))
      # compute rmse between predicted and actual after holdout procedure
      age.sum <- data.frame(Metric="RMSE", Dataset=exp$Dataset, nsub=length(unique(graphs$subjects)),
                            nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                            nroi=sqrt(dim(flat.gr$array)[2]), task="Age", thresh=exp$thr,
                            stat=rmse(age.res$true, age.res$pred), embed=embed, null=var(age.res$true))

      sex.res <- do.call(rbind, lapply(unique(graphs$subjects), function(sub) {
        training.set <- which(graphs$subjects != sub)  # hold out same-subjects from training set
        testing.set <- which(graphs$subjects == sub)  # validate over all scans for this subject
        trained.sex.rf <- randomForest(embed.graphs[training.set,], y=factor(pheno.scans$SEX[training.set]))
        preds.sex.rf <- predict(trained.sex.rf, embed.graphs[testing.set,])
        return(data.frame(true=as.numeric(as.character(pheno.scans$SEX[testing.set])),
                          pred=as.numeric(as.character(preds.sex.rf))))
      }))

      sex.sum <- data.frame(Metric="MR", Dataset=exp$Dataset, nsub=length(unique(graphs$subjects)),
                            nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                            nroi=sqrt(dim(flat.gr$array)[2]), task="Sex", thresh=exp$thr,
                            stat=mean(sex.res$true != sex.res$pred), embed=embed,
                            null=min(sapply(unique(pheno.scans$SEX), function(sex) mean(pheno.scans$SEX == sex))))

      return(rbind(age.sum, sex.sum))
    }, error=function(e) {return(NULL)})
  }))

  result <- list(statistics=res, problem=task.res)
  saveRDS(result, file.path(rf.res.path, paste0("rf_dset-", exp$Dataset, "_thr-", exp$thr*1000, ".rds")))

  return(result)
}, mc.cores=no_cores)

robj <- list(statistics=do.call(rbind, lapply(rf.results, function(r) r$statistics)),
             problem=do.call(rbind, lapply(rf.results, function(r) r$problem)))

saveRDS(robj, file.path(opath, "rf_fmri_results.rds"))


