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

fmri.path <- '/mnt/nfs2/MR/cpac_3-9-2/'
pheno.path <- '/mnt/nfs2/MR/all_mr/phenotypic/'
opath <- '.data/real/'


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
      vertices <- union(vertices, V(tgr))
    }
  }

  vertices <- order(vertices)
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

dsets <- list.dirs(path=fmri.path, recursive=FALSE)

# atlas_opts <- c("C", "D")
# names(atlas_opts) <- c("cc2", "des")
atlas_opts <- c("D")
names(atlas_opts) <- c("des")

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

  lapply(names(atlas_opts), function(atlas) {
    list(Dataset=dset.key, Parcellation=atlas_opts[atlas],
         dat.path=file.path(fmri.path, dset_name, "graphs", sprintf("FSL_nff_nsc_gsr_%s", atlas)),
         pheno.path=file.path(pheno.path, paste(dset.key, "_phenotypic_data.csv", sep="")),
         sub.pos = sub.pos)
  })
}))


trng <- seq(0, 1, .01)

stats <- list(discr.os, anova.os, icc.os, i2c2.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2")

#================
#
# fMRI Driver
#
#================

# range of thresholds to try
# multicore apply over dataset
mclapply(experiments, function(exp) {
  graphs <- cpac.open_graphs(exp$dat.path, rtype="array", dataset_id=exp$Dataset,
                             atlas_id=exp$Parcellation, sub_pos = exp$sub.pos, flatten = TRUE)
  res.thresh <- lapply(trng, function(thr) {
    test <- graphs
    # threshold the graphs
    test$graphs[test$graphs < thr] <- 0
    test$graphs[test$graphs >= thr] <- 1
    # run statistics
    res <- do.call(rbind, lapply(names(stats), function(stat) {
      tryCatch({
        return(data.frame(Dataset=exp$Dataset, thresh=thr, alg=stat,
                          nses=length(unique(test$sessions)), nscans=dim(test$graphs)[1],
                          nroi=sqrt(dim(test$graphs)[2]), nsub=length(unique(test$subjects)),
                          stat=do.call(stats[[stat]], list(test$graphs, test$subjects))))
      }, error=function(e) {return(NULL)})
    }))

    pheno.dat <- read.csv(exp$pheno.path)
    pheno.dat$AGE_AT_SCAN_1 <- as.numeric(as.character(pheno.dat$AGE_AT_SCAN_1))
    pheno.dat <- pheno.dat[!duplicated(pheno.dat$SUBID),]
    pheno.dat <- pheno.dat[, c("SUBID", "AGE_AT_SCAN_1", "SEX")]
    pheno.scans <- pheno.dat[sapply(as.numeric(test$subjects), function(x) which(x == pheno.dat$SUBID)),]
    if (exp$Dataset == "KKI2009") {
      pheno.scans$SEX <- as.factor((pheno.scans$SEX == "M") + 1)
    }
    problems <- lapply(1:10, function(r) {
      # aggregate results across all subjects; report RMSE at current r
      age.res <- do.call(rbind, lapply(unique(test$subjects), function(sub) {
        training.set <- which(test$subjects != sub)  # hold out same-subjects from training set
        testing.set <- which(test$subjects == sub)  # validate over all scans for this subject
        # predict for held-out subject
        age.preds <- knn.reg(test$graphs[training.set,], test$graphs[testing.set,],
                         as.numeric(pheno.scans$AGE_AT_SCAN_1[training.set]), k=r)
        return(data.frame(k=r, true=pheno.scans$AGE_AT_SCAN_1[testing.set],
                          pred=age.preds$pred))
      }))
      # compute rmse between predicted and actual after holdout procedure
      age.sum <- data.frame(k=r, Metric="RMSE", Dataset=exp$Dataset, nsub=length(unique(test$subjects)),
                            nses=length(unique(test$sessions)), nscans=dim(test$graphs)[1],
                            nroi=sqrt(dim(test$graphs)[2]), task="Age", thresh=thr,
                            stat=rmse(age.res$true, age.res$pred), null=NaN)

      sex.res <- do.call(rbind, lapply(unique(test$subjects), function(sub) {
        training.set <- which(test$subjects != sub)  # hold out same-subjects from training set
        testing.set <- which(test$subjects == sub)  # validate over all scans for this subject
        sex.preds <- knn(test$graphs[training.set,], test$graphs[testing.set,],
                         as.factor(pheno.scans$SEX[training.set]), k=r)
        return(data.frame(k=r, true=as.numeric(as.character(pheno.scans$SEX[testing.set])),
                          pred=as.numeric(as.character(sex.preds))))
      }))

      sex.sum <- data.frame(k=r, Metric="MR", Dataset=exp$Dataset, nsub=length(unique(test$subjects)),
                            nses=length(unique(test$sessions)), nscans=dim(test$graphs)[1],
                            nroi=sqrt(dim(test$graphs)[2]), task="Sex", thresh=thr,
                            stat=mean(sex.res$true != sex.res$pred),
                            null=mean(pheno.scans$SEX == 1))
      return(list(age=age.sum, sex=sex.sum))
    })
    age.agg <- do.call(rbind, lapply(problems, function(problem) problem$age))
    sex.agg <- do.call(rbind, lapply(problems, function(problem) problem$sex))

    return(list(statistics=res, problem=rbind(age.agg, sex.agg)))
  })

  robj <- list(statistics=do.call(rbind, lapply(res.thresh, function(r) r$statistics)),
               problem=do.call(rbind, lapply(res.thresh, function(r) r$problem)))

  saveRDS(robj, file.path(exp$dat.path, "knn_results.rds"))
  return(robj)
}, mc.cores=no_cores)

results <- list(statistics=do.call(rbind, lapply(fmri.results, function(res) res$statistics)),
                problem=do.call(rbind, lapply(fmri.results, function(res) res$problem)))
saveRDS(results, file.path(opath, "knn_fmri_results.rds"))


#############################################################
#############################################################
#############################################################

require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(igraph)
require(fmriutils)
require(reshape2)
require(stringr)
require(caret)

dmri.path <- '/data/nkienh/dmri/desikan/'
pheno.path <- '/data/all_mr/phenotypic/NKI1_phenotypic_data.csv'
opath <- '~/Documents/research/Rpackages/mgc/docs/discriminability/paper/data/real/'

no_cores = detectCores()/2

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

fmriu.io.open_graphs <- function(fnames, dataset_id="", atlas_id="",
                                 fmt='elist', verbose=FALSE, rtype='list', flatten=FALSE,
                                 rem.diag=TRUE) {
  if (! (fmt %in% c('adj', 'elist'))) {
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
      vertices <- union(vertices, V(tgr)$name)
    }
  }

  vertices <- vertices[ordered(as.numeric(vertices))]
  counter <- 1
  for (i in 1:length(fnames)) {
    basename <- basename(fnames[i])     # the base name of the file
    if (verbose) {
      print(paste('Loading', basename, '...'))
    }
    igr <- tryCatch({
      igraph::read_graph(fnames[i], format=fmt, predef=vertices) # read the graph from the filename, ordering by the vertices we found previously
    }, error = function(e) {
      return(NaN)
    })

    if (is.igraph(igr)) {
      tgr <- get.adjacency(igr, type="both", attr="weight", sparse=FALSE) # convert graph to adjacency matrix
      tgr[is.nan(tgr)] <- 0  # missing entries substituted with 0s
      if (rem.diag) {
        diag(tgr) <- 0
      }
      colnames(tgr) <- V(igr); rownames(tgr) <- V(igr)
      gr[[basename]] <-t(tgr)
      subjects[counter] <- str_extract(basename, 'sub(.?)+?(?=_)')
      sessions[counter] <- str_extract(basename, 'ses(.?)+?(?=_)')
      tasks[counter] <- str_extract(basename, 'task(.?)+?(?=_)')
      counter <- counter + 1
    }
  }

  dataset <- rep(dataset_id, counter - 1)
  atlas <- rep(atlas_id, counter - 1)
  subjects <- subjects[1:counter - 1]
  sessions <- sessions[1:counter - 1]
  tasks <- tasks[1:counter - 1]

  if (rtype == 'array') {
    aro <- fmriu.list2array(gr, flatten=flatten)
    gr <- aro$array
    dataset <- dataset[aro$incl_ar]
    atlas <- atlas[aro$incl_ar]
    subjects <- subjects[aro$incl_ar]
    sessions <- sessions[aro$incl_ar]
    tasks <- tasks[aro$incl_ar]
  }
  return(list(graphs=gr, dataset=dataset, atlas=atlas, subjects=subjects,
              sessions=sessions, tasks=tasks))
}


nthresh <- 100

pheno.dat <- read.csv(pheno.path)

stats <- list(discr.os, anova.os, icc.os, i2c2.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2")

#================
#
# dMRI Driver
#
#================
# open the graphs
graphs <- fmriu.io.open_graphs(dmri.path, rtype="array", dataset_id="NKI24", atlas_id="Desikan", flatten = TRUE)
# qet the 0 -> 1 quantiles in .01 increments
trng <- quantile(graphs$graphs[graphs$graphs != 0], seq(0, 1, .01))
pheno.linked <- pheno.dat
pheno.linked <- pheno.linked[!duplicated(pheno.dat$SUBID),]
# Class 1 if they have a current drug addiction or mental illness diagnosis; 0 otherwise
pheno.linked$Lifestyle <- as.numeric(pheno.linked$LIFETIME_DX_1 != "#" | pheno.linked$LIFETIME_DX_2 != "#" |
                                       pheno.linked$LIFETIME_DX_3 != "#" | pheno.linked$LIFETIME_DX_4 != "#" |
                                       pheno.linked$LIFETIME_DX_5 != "#")
pheno.linked$Mental <- pheno.linked$CURRENT_DIAGNOSIS
pheno.linked$Age <- pheno.linked$AGE_AT_SCAN_1

# mclapply over it
dmri.results <- mclapply(trng, function(thr) {
  test <- graphs
  # threshold the graphs
  test$graphs[test$graphs < thr] <- 0
  test$graphs[test$graphs >= thr] <- 1
  # run statistics
  res <- do.call(rbind, lapply(names(stats), function(stat) {
    tryCatch({
      return(data.frame(thresh=thr, alg=stat, stat=do.call(stats[[stat]], list(test$graphs, test$subjects))))
    }, error=function(e) {return(NULL)})
  }))
  test$subjects <- gsub(".*-", "", test$subjects)
  pheno.dat <- pheno.linked[sapply(as.numeric(test$subjects), function(x) which(x == pheno.linked$SUBID)),]

  fit.age <- train(V1~.,
                   method="knn",
                   tuneGrid=expand.grid(k=1:10),
                   trControl=trainControl(method  = "LOOCV"),
                   metric="RMSE",
                   data=data.frame(cbind(V1=as.numeric(pheno.dat$Age), test$graphs))
  )
  dat <- data.frame(cbind(V1=pheno.dat$Mental, test$graphs))
  dat$V1 <- as.factor(make.names(as.factor(dat$V1)))
  fit.mental <- train(V1~.,
                      method="knn",
                      tuneGrid=expand.grid(k=1:10),
                      trControl=trainControl(method="LOOCV", classProbs=TRUE),
                      metric="Accuracy",
                      data=dat
  )
  dat <- data.frame(cbind(V1=pheno.dat$Lifestyle, test$graphs))
  dat$V1 <- as.factor(make.names(as.factor(dat$V1)))
  fit.lifestyle <- train(V1 ~ .,
                         method="knn",
                         tuneGrid=expand.grid(k=1:10),
                         trControl=trainControl(method="LOOCV", classProbs=TRUE),
                         metric="Accuracy",
                         data=dat
  )
  fit.age$results <- melt(fit.age$results, id="k")
  names(fit.age$results) <- c("k", "Metric", "value"); fit.age$results$task = "Age"

  fit.mental$results <- fit.mental$results[, c("k", "Accuracy")]
  fit.mental$results <- melt(fit.mental$results, id="k")
  names(fit.mental$results) <- c("k", "Metric", "value"); fit.mental$results$task = "Mental"

  fit.lifestyle$results <- fit.lifestyle$results[, c("k", "Accuracy")]
  fit.lifestyle$results <- melt(fit.lifestyle$results, id="k")
  names(fit.lifestyle$results) <- c("k", "Metric", "value"); fit.lifestyle$results$task = "Lifestyle"

  results <- rbind(fit.age$results, fit.mental$results, fit.lifestyle$results)
  results$thresh <- thr
  return(list(statistics=res, problem=results))
}, mc.cores=no_cores/2)

results <- list(statistics=do.call(rbind, lapply(dmri.results, function(res) res$statistics)),
                problem=do.call(rbind, lapply(dmri.results, function(res) res$problem)))
saveRDS(results, file.path(opath, "knn_dmri_results.rds"))

