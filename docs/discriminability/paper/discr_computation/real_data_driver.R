# ==============================
# fMRI Driver
# ==============================

require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(igraph)
require(stringr)
require(fmriutils)

no_cores = detectCores() - 31

ipath <- "./"
opath <- "./"

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

dsets <- list.dirs(path=ipath, recursive=FALSE)

dsets <- dsets[!(dsets %in% c(".//MPG1", ".//BNU3"))]
dsets <- dsets[dsets %in% c(".//NKI24_mx1400", ".//NKI24_mx645", ".//NKI24_std2500", ".//KKI2009")]
experiments <- do.call(c, lapply(dsets, function(dset) {
  dset_name = basename(dset)
  do.call(c, lapply(names(reg_opts), function(reg) {
    do.call(c, lapply(names(freq_opts), function(freq) {
      do.call(c, lapply(names(scrub_opts), function(scrub) {
        do.call(c, lapply(names(gsr_opts), function(gsr) {
          lapply(names(atlas_opts), function(atlas) {
            list(Dataset=dset_name, Reg=reg_opts[reg],
                 FF=freq_opts[freq], Scr=scrub_opts[scrub], GSR=gsr_opts[gsr],
                 Parcellation=atlas_opts[atlas],
                 path=file.path(ipath, dset_name, "graphs", paste(reg, freq,scrub, gsr, atlas, sep="_")))
          })
        }))
      }))
    }))
  }))
}))

results <- mclapply(experiments, function(exp) {
  if (grepl("NKI", exp$path)) {
    sub.pos <- 3
  } else if (grepl("KKI", exp$path)) {
    sub.pos <- 1
  } else {
    sub.pos <- 2
  }
  if (FALSE) {#file.exists(file.path(exp$path, "discr_results.rds"))) {
    res <- readRDS(file.path(exp$path, "discr_results.rds"))
  } else {
    graphs <- cpac.open_graphs(exp$path, rtype="array", dataset_id=exp$Dataset, atlas_id=exp$Parcellation,
                                     sub_pos=sub.pos, flatten = TRUE)
    res <- lapply(names(graph.xfm), function(xfm) {
      tryCatch({
        test <- graphs
        test$graphs <- t(do.call(apply, list(test$graphs, 1, graph.xfm[[xfm]], min(test$graphs[test$graphs != 0]))))
        D=discr.distance(test$graphs)
        res <- discr.stat(test$graphs, test$subjects)
        return(list(D=D, data=data.frame(Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
                          Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation,
                          xfm=xfm, nsub=length(unique(test$subjects)),
                          nses=length(unique(test$sessions)), nscans=dim(test$graphs)[1], nroi=sqrt(dim(test$graphs)[2]),
                          discr=res$discr), subjects=graphs$subjects))
      }, error=function(e) {
        list(D=D, data=data.frame(Dataset=exp$Dataset, Reg=exp$Reg, FF=exp$FF,
                   Scr=exp$Scr, GSR=exp$GSR, Parcellation=exp$Parcellation,
                   xfm=xfm, nsub=length(unique(test$subjects)),
                   nses=length(unique(test$sessions)), nscans=dim(test$graphs)[1], nroi=sqrt(dim(test$graphs)[2]),
                   discr=NaN), subjects=graphs$subjects)
        })
    })
    saveRDS(res, file.path(exp$path, "discr_results.rds"))
  }
  return(do.call(rbind, lapply(res, function(x) x$data)))
  }, mc.cores=no_cores)

results.fmri <- do.call(rbind, results)
saveRDS(results.fmri, file.path(opath, paste('discr_fmri_results', '.rds', sep="")))

# =====================
# dMRI Driver
# =====================

require(lolR)
require(MASS)
library(parallel)
require(mgc)
require(ICC)
require(I2C2)
require(igraph)
require(stringr)
require(fmriutils)
no_cores = detectCores() - 1

ipath <- "./"
opath <- "./"


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
  return(log(x + min.x/exp(1)))
}

graph.xfm <- list(nofn, ptr, log.xfm)
names(graph.xfm) <- c("N", "P", "L")


dsets <- list.dirs(path=ipath, recursive=FALSE)

dsets <- dsets[dsets %in% c(".//BNU1", ".//HNU1", ".//SWU4", ".//KKI2009", ".//NKI24")]
experiments <- do.call(c, lapply(dsets, function(dset) {
  dset_name <- basename(dset)
  at_dirs <- list.dirs(file.path(ipath, dset_name, "ndmg_0-0-48", "graphs"))
  lapply(at_dirs, function(atlas) {
    if (basename(atlas) != "graphs") {
      if (length(list.files(atlas, pattern="\\.ssv"))) {
        return(list(Dataset=dset_name, Parcellation=basename(atlas), path=atlas))
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
}))
experiments[sapply(experiments, is.null)] <- NULL

results <- mclapply(experiments, function(exp) {
  graphs <- fmriu.io.open_graphs(exp$path, rtype="array", dataset_id = exp$Dataset, atlas_id = exp$Parcellation,
                                 flatten = TRUE, rem.diag=TRUE)
  res <- lapply(names(graph.xfm), function(xfm) {
    tryCatch({
      test <- graphs
      test$graphs <- t(do.call(apply, list(test$graphs, 1, graph.xfm[[xfm]], min(test$graphs[test$graphs != 0]))))
      D <- discr.distance(test$graphs)
      res <- discr.stat(D, test$subjects)
      return(list(D=D, data=data.frame(Dataset=exp$Dataset, Parcellation=exp$Parcellation,
                       xfm=xfm,
                       nsub=length(unique(graphs$subjects)), nses=length(unique(graphs$sessions)),
                       nscans=dim(graphs$graphs)[1], nroi=sqrt(dim(graphs$graphs))[2], discr=res$discr),
                  subjects=graphs$subjects))
    }
    , error=function(e) {print(e); return(list(D=D, data=data.frame(Dataset=exp$Dataset, Parcellation=exp$Parcellation,
                                                     xfm=xfm,
                                                     nsub=length(unique(graphs$subjects)), nses=length(unique(graphs$sessions)),
                                                     nscans=dim(graphs$graphs)[1], nroi=sqrt(dim(graphs$graphs))[2], discr=res$discr),
                                               subjects=graphs$subjects)
    )})
    })
  saveRDS(res, file.path(exp$path, "discr_results.rds"))
  return(do.call(rbind, lapply(res, function(x) x$data)))
  }, mc.cores=no_cores)

results.dmri <- do.call(rbind, results)
saveRDS(results.dmri, file.path(opath, paste('discr_dmri_results', '.rds', sep="")))
