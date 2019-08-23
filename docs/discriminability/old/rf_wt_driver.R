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
# load/source MASE code
mase.path <- './mase/R/'
mase.files <- list.files(mase.path)
mase.files <- mase.files[mase.files %in% c("mase.R", "omnibus-embedding.R", "getElbows.R")]
sapply(mase.files, function(x) source(file.path(mase.path, x)))

fmri.path <- '/mnt/nfs2/MR/cpac_3-9-2/'
pheno.path <- '/mnt/nfs2/MR/all_mr/phenotypic/'
#fmri.path <- '/cis/project/ndmg/eric/discriminability/cpac_3-9-2/'
#pheno.path <- '/cis/project/ndmg/eric/discriminability/phenotypic/'
#fmri.path <- '/data/cpac_3-9-2/'
#pheno.path <- '/data/all_mr/phenotypic/'
opath <- './data/real/'
no_cores <- parallel::detectCores() - 10

mgc.testt <- function(x, y, R=1000) {
  result <- mgc.test(X=x, Y=y, rep=R)
  return(list(p.value=result$pMGC, statistic=result$statMGC))
}

# dependence test methods
dep.tests <- list(mgc=mgc.testt, dcor=dcor.test)


# dependence test methods
dep.tests <- list(mgc=mgc.testt, dcor=dcor.test)

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

nofn <- function(x, ...) {
  return(x)
}

ptr <- function(x, ...) {
  nz <- x[x != 0]
  r <- rank(nz)*2/(length(nz) + 1)
  x[x != 0] <- r
  x <- (x - min(x))/(max(x) - min(x))
  return(x)
}

log.xfm <- function(x, min.x) {
  return(log(x + min.x/exp(2)))
}

reg_opts <- c("A", "F")
names(reg_opts) <- c("ANT", "FSL")
freq_opts <- c("F", "N")
names(freq_opts) <- c("frf", "nff")
scrub_opts <- c("S", "N")
names(scrub_opts) <- c("scr", "nsc")
gsr_opts <- c("G", "N")
names(gsr_opts) <- c("gsr", "ngs")
atlas_opts <- c("A", "C", "D", "H")
names(atlas_opts) <- c("aal", "cc2", "des", "hox")

graph.xfm <- list(nofn, ptr, log.xfm)
names(graph.xfm) <- c("N", "P", "L")

dsets <- list.dirs(path=fmri.path, recursive=FALSE)

dsets <- dsets[!(dsets %in% c(".//MPG1", ".//BNU3"))]
#dsets <- dsets[dsets %in% c(".//NKI24_mx1400", ".//NKI24_mx645", ".//NKI24_std2500", ".//KKI2009")]
experiments <- do.call(c, lapply(dsets, function(dset) {
  dset_name = basename(dset)
  do.call(c, lapply(names(reg_opts), function(reg) {
    do.call(c, lapply(names(freq_opts), function(freq) {
      do.call(c, lapply(names(scrub_opts), function(scrub) {
        do.call(c, lapply(names(gsr_opts), function(gsr) {
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
            list(Dataset=dset.key, Reg=reg_opts[reg],
                 FF=freq_opts[freq], Scr=scrub_opts[scrub], GSR=gsr_opts[gsr], Parcellation=atlas_opts[atlas],
                 dat.path=file.path(fmri.path, dset_name, "graphs", paste(reg, freq,scrub, gsr, atlas, sep="_")),
                 pheno.path=file.path(pheno.path, paste(dset.key, "_phenotypic_data.csv", sep="")),
                 sub.pos = sub.pos)
          })
        }))
      }))
    }))
  }))
}))


stats <- list(discr.os, anova.os, icc.os, i2c2.os, manova.os)
names(stats) <- c("Discr", "ANOVA", "ICC", "I2C2", "MANOVA")

graph.xfms <- list(nofn, ptr, log.xfm)
names(graph.xfms) <- c("N", "P", "L")

dep.res.path <- file.path(opath, 'dep_results_wt')
dir.create(opath)
dir.create(dep.res.path)

#================
#
# fMRI Driver
#
#================

# range of thresholds to try
# multicore apply over dataset
dep.results <- mclapply(experiments, function(exp) {
  o.path <- file.path(dep.res.path, paste0("dep_dset-", exp$Dataset, "_",
                                          paste0(exp$Reg, exp$FF, exp$Scr, exp$GSR, exp$Parcellation), ".rds"))
  tryCatch({
    #if (file.exists(o.path)) {
    #  return(readRDS(o.path))
    #} else {
      graphs <- cpac.open_graphs(exp$dat.path, dataset_id=exp$Dataset,
                                 atlas_id=exp$Parcellation, sub_pos = exp$sub.pos, flatten=FALSE)

      print(sprintf("Dataset: %s, Reg=%s, FF=%s, Scr=%s, GSR=%s, Parc=%s", exp$Dataset, Reg=exp$Reg, FF=exp$FF,
                    Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation))

      # flatten for statistics
      result <- lapply(names(graph.xfms), function(graph.xfm) {
        test <- graphs
        min.gr <- min(sapply(test$graphs, function(gr) min(gr[gr != 0])))
        test$graphs <- lapply(test$graphs, function(gr) do.call(graph.xfms[[graph.xfm]], list(gr, min.gr)))
        flat.gr <- fmriu.list2array(test$graphs, flatten=TRUE)

        stat.res <- do.call(rbind, lapply(names(stats), function(stat) {
          tryCatch({
            return(data.frame(Dataset=exp$Dataset, alg=stat, Reg=exp$Reg, FF=exp$FF,
                              Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation,
                              xfm=graph.xfm, nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                              nroi=sqrt(dim(flat.gr$array)[2]), nsub=length(unique(graphs$subjects)),
                              stat=do.call(stats[[stat]], list(flat.gr$array, graphs$subjects))))
          }, error=function(e) {return(NULL)})
        }))

        graphs.embedded <- list(
          raw=t(simplify2array(lapply(test$graphs, function(x) as.vector(x)))),
          mase=t(simplify2array(lapply(mase(test$graphs)$R, function(x) as.vector(x))))
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

        # task.res <- do.call(rbind, lapply(names(graphs.embedded), function(embed) {
        #   tryCatch({
        #     embed.graphs <- graphs.embedded[[embed]]
        #     # aggregate results across all subjects; report RMSE at current r
        #     age.res <- do.call(rbind, lapply(unique(graphs$subjects), function(sub) {
        #       training.set <- which(graphs$subjects != sub)  # hold out same-subjects from training set
        #       testing.set <- which(graphs$subjects == sub)  # validate over all scans for this subject
        #       # predict for held-out subject
        #       trained.age.rf <- randomForest(embed.graphs[training.set,], y=as.numeric(pheno.scans$AGE_AT_SCAN_1[training.set]))
        #       preds.age.rf <- predict(trained.age.rf, embed.graphs[testing.set,])
        #       return(data.frame(true=pheno.scans$AGE_AT_SCAN_1[testing.set],
        #                         pred=preds.age.rf, subject=sub))
        #     }))
        #     # compute rmse between predicted and actual after holdout procedure
        #     age.sum <- data.frame(Metric="RMSE", Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
        #                           Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation, xfm=graph.xfm,
        #                           nsub=length(unique(graphs$subjects)),
        #                           nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
        #                           nroi=sqrt(dim(flat.gr$array)[2]), task="Age",
        #                           stat=rmse(age.res$true, age.res$pred), embed=embed, null=var(age.res$true))
        #
        #     sex.res <- do.call(rbind, lapply(unique(graphs$subjects), function(sub) {
        #       training.set <- which(graphs$subjects != sub)  # hold out same-subjects from training set
        #       testing.set <- which(graphs$subjects == sub)  # validate over all scans for this subject
        #       trained.sex.rf <- randomForest(embed.graphs[training.set,], y=factor(pheno.scans$SEX[training.set]))
        #       preds.sex.rf <- predict(trained.sex.rf, embed.graphs[testing.set,])
        #       return(data.frame(true=as.numeric(as.character(pheno.scans$SEX[testing.set])),
        #                         pred=as.numeric(as.character(preds.sex.rf))))
        #     }))
        #
        #     sex.sum <- data.frame(Metric="MR", Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
        #                           Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation, xfm=graph.xfm,
        #                           nsub=length(unique(graphs$subjects)),
        #                           nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
        #                           nroi=sqrt(dim(flat.gr$array)[2]), task="Sex",
        #                           stat=mean(sex.res$true != sex.res$pred), embed=embed,
        #                           null=min(sapply(unique(pheno.scans$SEX), function(sex) mean(pheno.scans$SEX == sex))))
        #
        #     return(rbind(age.sum, sex.sum))
        #   }, error=function(e) {return(NULL)})
        # }))


        dcor.res <- do.call(rbind, lapply(names(graphs.embedded), function(embed) {
          embed.graphs <- graphs.embedded[[embed]]
          # aggregate results across all subjects; report RMSE at current r
          return(do.call(rbind, lapply(names(dep.tests), function(dep) {
            tryCatch({
              Y.age <- as.numeric(pheno.scans$AGE_AT_SCAN_1)
              valid.idx.age <- (!is.na(Y.age) & !is.null(Y.age) & !is.nan(Y.age))
              Y.sex <- as.numeric(pheno.scans$SEX)
              valid.idx.sex <- (!is.na(Y.sex) & !is.null(Y.sex) & !is.nan(Y.sex))
              dep.age <- do.call(dep.tests[[dep]], list(x=embed.graphs[valid.idx.age,], y=Y.age[valid.idx.age], R=1000))
              dep.sex <- do.call(dep.tests[[dep]], list(x=embed.graphs[valid.idx.sex,], y=Y.sex[valid.idx.sex], R=1000))
              return(rbind(data.frame(Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
                                      Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation, xfm=graph.xfm,
                                      nsub=length(unique(graphs$subjects)),
                                      nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                                      nroi=sqrt(dim(flat.gr$array)[2]), task="Age",
                                      stat=dep.age$statistic, pval=dep.age$p.value, method=dep),
                           data.frame(Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
                                      Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation, xfm=graph.xfm,
                                      nsub=length(unique(graphs$subjects)),
                                      nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                                      nroi=sqrt(dim(flat.gr$array)[2]), task="Sex",
                                      stat=dep.sex$statistic, pval=dep.sex$p.value, method=dep)))
            }, error=function(e) {return(e)})
          })))
        }))

        return(list(statistics=stat.res, problem=task.res, dcor=dcor.res))
      })
      stat.result <- do.call(rbind, lapply(result, function(res) res$statistics))
      #task.result <- do.call(rbind, lapply(result, function(res) res$problem))
      dcor.result <- do.call(rbind, lapply(result, function(res) res$dcor))

      #result <- list(statistics=stat.result, problem=task.result, dcor=dcor.result)
      result <- list(statistics=stat.result, dcor=dcor.result)

      saveRDS(result, o.path)
      return(result)
    #}
  }, error=function(e) {return(NULL)})
}, mc.cores=no_cores)

robj <- list(statistics=do.call(rbind, lapply(dep.results, function(res) res$statistics)),
             #problem=do.call(rbind, lapply(dep.results, function(res) res$problem)),
             dcor=do.call(rbind, lapply(dep.results, function(res) res$dcor)))

saveRDS(robj, file.path(opath, "dep_wt_fmri_results.rds"))


