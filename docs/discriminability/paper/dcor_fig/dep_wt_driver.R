require(parallel)
require(mgc)
require(lolR)
require(I2C2)
require(ICC)
require(igraph)
require(reshape2)
require(stringr)
require(FNN)
require(Metrics)
require(randomForest)
require(rARPACK)
require(energy)
source('../simulations/shared_scripts.R')


fmri.path <- '/mnt/nfs2/MR/cpac_3-9-2/' # '/mnt/nfs2/MR/cpac_3-9-2/'
pheno.path <- '/mnt/nfs2/MR/all_mr/phenotypic/' # '/mnt/nfs2/MR/all_mr/phenotypic/'
#fmri.path <- '/cis/project/ndmg/eric/discriminability/cpac_3-9-2/'
#pheno.path <- '/cis/project/ndmg/eric/discriminability/phenotypic/'
#fmri.path <- '/data/cpac_3-9-2/'
#pheno.path <- '/data/all_mr/phenotypic/'
opath <- '../data/real/'
no_cores <- parallel::detectCores() - 1

dep.tests <- list(DCorr=ksamp.test)

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

nofn.xfm <- function(x, ...) {
  return(x)
}

ptr.xfm <- function(x, ...) {
  nz <- x[x != 0]
  r <- rank(nz)*2/(length(nz) + 1)
  x[x != 0] <- r
  x <- (x - min(x))/(max(x) - min(x))
  return(x)
}

log.xfm <- function(x, min.x=10^(-6)) {
  return(log(x + min.x/exp(2)))
}

stats <- list(Discr=discr.os, PICC=icc.os, I2C2=i2c2.os, Kernel=ksamp.os, FPI=fpi.os,
              DISCO=disco.os, HSIC=hsic.os)#, manova.os)

graph.xfms <- list(nofn.xfm, ptr.xfm, log.xfm)
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
dep.results <- mclapply(experiments, function(exper) {
  t0 = Sys.time()
  o.path <- file.path(dep.res.path, paste0("dep_dset-", exper$Dataset, "_",
                                          paste0(exper$Reg, exper$FF, exper$Scr, exper$GSR, exper$Parcellation), ".rds"))
  tryCatch({
    #if (file.exists(o.path)) {
    # return(readRDS(o.path))
    #} else if (file.exists(exper$pheno.path)) {
      graphs <- cpac.open_graphs(exper$dat.path, dataset_id=exper$Dataset,
                                 atlas_id=exper$Parcellation, sub_pos = exper$sub.pos, flatten=FALSE)

      print(sprintf("Dataset: %s, Reg=%s, FF=%s, Scr=%s, GSR=%s, Parc=%s", exper$Dataset, Reg=exper$Reg, FF=exper$FF,
                    Scr=exper$Scr, GSR=exper$GSR, Parcellation=exper$Parcellation))

      # flatten for statistics
      result <- lapply(names(graph.xfms), function(graph.xfm) {
        print(graph.xfm)
        test <- graphs
        min.gr <- min(sapply(test$graphs, function(gr) min(gr[gr != 0])))
        test$graphs <- lapply(test$graphs, function(gr) do.call(graph.xfms[[graph.xfm]], list(gr, min.x=min.gr)))
        flat.gr <- fmriu.list2array(test$graphs, flatten=TRUE)

        stat.res <- do.call(rbind, lapply(names(stats), function(stat) {
          tryCatch({
            return(data.frame(Dataset=exper$Dataset, alg=stat, Reg=exper$Reg, FF=exper$FF,
                              Scr=exper$Scr, GSR=exper$GSR, Parcellation=exper$Parcellation,
                              xfm=graph.xfm, nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                              nroi=sqrt(dim(flat.gr$array)[2]), nsub=length(unique(graphs$subjects)),
                              stat=do.call(stats[[stat]], list(X=flat.gr$array, Y=graphs$subjects,
                                                               Z=graphs$sessions, is.dist=FALSE))))
          }, error=function(e) {return(NULL)})
        }))

        pheno.dat <- read.csv(exper$pheno.path)
        pheno.dat$AGE_AT_SCAN_1 <- as.numeric(as.character(pheno.dat$AGE_AT_SCAN_1))
        pheno.dat <- pheno.dat[!duplicated(pheno.dat$SUBID),]
        pheno.dat <- pheno.dat[, c("SUBID", "AGE_AT_SCAN_1", "SEX")]
        matched.idx <- lapply(as.numeric(graphs$subjects), function(x) {
          match <- which(x == pheno.dat$SUBID)
          if (length(match) > 0) {
            return(match)
          } else {
            return(NULL)
          }
        })
        not.retain.idx <- sapply(matched.idx, is.null)

        matched.idx <- unlist(matched.idx[!not.retain.idx])
        test$graphs <- test$graphs[!not.retain.idx]
        pheno.scans <- pheno.dat[matched.idx,]
        if (exper$Dataset == "KKI2009") {
          pheno.scans$SEX <- as.factor((pheno.scans$SEX == "M") + 1)
        }

        graphs.embedded <- list(
          Raw=t(simplify2array(lapply(test$graphs, function(x) as.vector(as.matrix(x)))))#,
          # mase=t(simplify2array(lapply(mase(test$graphs)$R, function(x) as.vector(x))))
        )
        #graphs.embedded$dist <- g.ase(as.matrix(dist(graphs.embedded$mase)))$X

        dcor.res <- do.call(rbind, lapply(names(graphs.embedded), function(embed) {
          embed.graphs <- graphs.embedded[[embed]]
          # aggregate results across all subjects; report RMSE at current r
          return(do.call(rbind, lapply(names(dep.tests), function(dep) {
            tryCatch({
              Y.age <- as.numeric(pheno.scans$AGE_AT_SCAN_1)
              valid.idx.age <- (!is.na(Y.age) & !is.null(Y.age) & !is.nan(Y.age))
              Y.sex <- as.numeric(pheno.scans$SEX)
              valid.idx.sex <- (!is.na(Y.sex) & !is.null(Y.sex) & !is.nan(Y.sex))
              dep.age <- do.call(dep.tests[[dep]], list(embed.graphs[valid.idx.age,], Y.age[valid.idx.age], is.dist=FALSE, nrep=1000L))
              dep.sex <- do.call(dep.tests[[dep]], list(embed.graphs[valid.idx.sex,], Y.sex[valid.idx.sex], is.dist=FALSE, nrep=1000L))
              return(rbind(data.frame(Dataset=exper$Dataset, Reg=exper$Reg, FF=exper$FF,
                                      Scr=exper$Scr, GSR=exper$GSR, Parcellation=exper$Parcellation, xfm=graph.xfm,
                                      nsub=length(unique(graphs$subjects)),
                                      nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                                      nroi=sqrt(dim(flat.gr$array)[2]), task="Age", embed=embed,
                                      stat=dep.age$statistic, method=dep, pvalue=dep.age$pvalue),
                           data.frame(Dataset=exper$Dataset, Reg=exper$Reg, FF=exper$FF,
                                      Scr=exper$Scr, GSR=exper$GSR, Parcellation=exper$Parcellation, xfm=graph.xfm,
                                      nsub=length(unique(graphs$subjects)),
                                      nses=length(unique(graphs$sessions)), nscans=dim(flat.gr$array)[1],
                                      nroi=sqrt(dim(flat.gr$array)[2]), task="Sex", embed=embed,
                                      stat=dep.sex$statistic, method=dep, pvalue=dep.sex$pvalue)))
            }, error=function(e) {print(e)})
          })))
        }))

        return(list(statistics=stat.res, dcor=dcor.res))
      })
      stat.result <- do.call(rbind, lapply(result, function(res) res$statistics))
      dcor.result <- do.call(rbind, lapply(result, function(res) res$dcor))

      result <- list(statistics=stat.result, dcor=dcor.result)

      saveRDS(result, o.path)
      return(result)
    #}
  }, error=function(e) {
    print(sprintf("Dataset: %s, Reg=%s, FF=%s, Scr=%s, GSR=%s, Parc=%s, ERR=%s", exper$Dataset, Reg=exper$Reg, FF=exper$FF,
                  Scr=exper$Scr, GSR=exper$GSR, Parcellation=exper$Parcellation, e))
  })
  t1 = Sys.time()
  print(t1-t0)
}, mc.cores=no_cores, mc.preschedule=FALSE)

dep.results.retain.idx <- sapply(dep.results, function(res) {
  tryCatch({
  if (ncol(res$statistics) != 13 || ncol(res$dcor) != 16) {
    return(FALSE)
  } else {
    return(TRUE)
  }}, error=function(e) {return(FALSE)})
})
dep.results.purged <- dep.results[dep.results.retain.idx]

robj <- list(statistics=do.call(rbind, lapply(dep.results.purged, function(res) {
  tryCatch({res$statistics}, error=function(e) {return(NULL)})
  })),
  dcor=do.call(rbind, lapply(dep.results.purged, function(res)  {
    tryCatch({res$dcor}, error=function(e) {return(NULL)})
  })))

saveRDS(list(processed=robj, raw=dep.results), file.path(opath, "dep_wt_fmri_results.rds"))


