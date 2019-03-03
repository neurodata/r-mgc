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

fmri.path <- '/mnt/nfs2/MR/cpac_3-9-2/'
opath <- './data/real/'

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

dsets <- dsets[!sapply(dsets, function(dset) {grepl("BNU3", dset) | grepl("MPG1", dset)})]
experiments <- lapply(dsets, function(dset) {
  dset_name = basename(dset)
  do.call(c, lapply(names(reg_opts), function(reg) {
    do.call(c, lapply(names(freq_opts), function(freq) {
      do.call(c, lapply(names(scrub_opts), function(scrub) {
        do.call(c, lapply(names(gsr_opts), function(gsr) {
          do.call(c, lapply(names(atlas_opts), function(atlas) {
            dat=readRDS(file.path(fmri.path, dset_name, "graphs", paste(reg, freq,scrub, gsr, atlas, sep="_"), "discr_results.rds"))
            lapply(dat, function(d) {
              list(D=d$D, data=d$dat, subjects=d$subjects, opath=file.path(fmri.path, dset_name, "graphs", "pval_results.rds"))
            })
          }))
        }))
      }))
    }))
  }))
})

fmri.results <- lapply(experiments, function(dset.exp) {
  print(print(sprintf("%s Dataset...", dset.exp[[1]]$data$Dataset)))
  exp.sets <- do.call(c, lapply(1:length(dset.exp), function(i) {
    lapply(1:length(dset.exp), function(j) {
      return(list(i=i,j=j))
    })
  }))
  exp.res <- mclapply(1:length(exp.sets), function(idx) {
    ex=exp.sets[[idx]]
    if (idx %% 1000 == 0) {
      print(sprintf("Comparison %d of %d for %.2f percent...", idx,
                    length(exp.sets), 100*idx/length(exp.sets)))
    }
    pval <- discr.test.two_sample(dset.exp[[ex$i]]$D, dset.exp[[ex$j]]$D, 
                                  ids=dset.exp[[ex$i]]$subjects, nperm=1000)
    return(list(i=ex$i, j=ex$j, pval=pval$pval))
  }, mc.cores=no_cores)
  
  pval.mtx <- array(NaN, dim=c(length(dset.exp), length(dset.exp)))
  
  for (res in exp.res) {
    pval.mtx[res$i, res$j] <- res$pval
  }
  
  pipes <- do.call(rbind, lapply(dset.exp, function(dset.pipe) {
    return(dset.pipe$data)
  }))
  pipes$Names <- as.character(apply(pipes[, c("Reg", "FF", "Scr", "GSR", "Parcellation", "xfm")], 
                                    c(1), function(x) {paste(x, collapse="")}))
  rownames(pval.mtx) <- pipes$Names
  colnames(pval.mtx) <- pipes$Names
  robj <- list(pvals=pval.mtx, pipes=pipes)
  saveRDS(robj, dset.exp[[1]]$opath)
  return(robj)
})
names(fmri.results) <- dsets
saveRDS(fmri.results, file.path(opath, "pval_results.rds"))

