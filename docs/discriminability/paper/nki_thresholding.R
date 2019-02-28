require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(igraph)
require(fmriutils)
require(reshape2)
require(caret)
require(stringr)

fmri.path <- '/data/nkienh/fmri/desikan/'
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


nthresh <- 100

pheno.dat <- read.csv(pheno.path)

stats <- list(discr.os, anova.os, icc.os, i2c2.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2")

#================
#
# fMRI Driver
#
#================
# open the graphs
graphs <- cpac.open_graphs(fmri.path, rtype="array", dataset_id="NKI24", atlas_id="Desikan", sub_pos = 3, flatten = TRUE)
# qet the 0 -> 1 quantiles in .01 increments
trng <- seq(0, 1, .01)
pheno.linked <- pheno.dat
pheno.linked <- pheno.linked[!duplicated(pheno.dat$SUBID),]
# Class 1 if they have a current drug addiction or mental illness diagnosis; 0 otherwise
pheno.linked$Lifestyle <- as.numeric(pheno.linked$LIFETIME_DX_1 != "#" | pheno.linked$LIFETIME_DX_2 != "#" |
                                       pheno.linked$LIFETIME_DX_3 != "#" | pheno.linked$LIFETIME_DX_4 != "#" |
                                       pheno.linked$LIFETIME_DX_5 != "#")
pheno.linked$Mental <- pheno.linked$CURRENT_DIAGNOSIS
pheno.linked$Age <- pheno.linked$AGE_AT_SCAN_1

# mclapply over it
fmri.results <- mclapply(trng, function(thr) {
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
  print(res)

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

