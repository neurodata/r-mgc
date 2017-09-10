#setwd("/users/ylee/multitime")

#source("/users/ylee/newMGC/MGCDistTransform.R")
#source("/users/ylee/newMGC/FindLargestRectangles.R")
#source("/users/ylee/newMGC/MGCLocalCorr.R")
#source("/users/ylee/newMGC/MGCSampleStat.R")
#source("/users/ylee/newMGC/MGCPermutationTest.R")

source("Auxiliary/connect_graph.R")
source("Auxiliary/notsimple3block.R")
source("Auxiliary/diffmaps.R")
source("Auxiliary/dmaps.R")
source("Auxiliary/elbowmap.R")
source("Auxiliary/getElbows.R")
source("Auxiliary/diffusionMGC.R")


library(igraph)
library(network)
library(MASS)
library(xtable)
library(parallel)
library(foreach)
library(HHG)
library(lattice)
library(DTMCPack)
library(Matrix)
library(energy)
library(Rlab)
library(VGAM)
library(amen)
library(doParallel)
library(ecodist)
library(SDMTools)
library(MGC)

######## setting ################
popn = 100
p = 0.5; q = 0.2; r = 0.3
n.perm = 200; n.iter = n.perm
nr = 200 # the number of independent iterations

mgc.result = list(); mcorr.result = list(); hhg.result = list()

for(i in 1:nr){
  print(i)
		set.seed(i)

    G = notsimple3block(popn, p, q, r)
  	X = V(G)$outcome

  	A = as.matrix(get.adjacency(G))

  	P = A  / rowSums(A)
    diffusion.q  =  min( max(getElbows(print.lambda(P, times = 3)[[1]], plot = FALSE, n = 3, threshold = 0)), popn-1)

    mgc.result[[i]] =  NetworkTest.diffusion.stat(G, X, option = 1, diffusion = TRUE, t.range = 3, n.perm = n.perm, q = diffusion.q)
    mcorr.result[[i]] =  NetworkTest.diffusion.stat(G, X, option = 2, diffusion = TRUE, t.range = 3, n.perm = n.perm, q = diffusion.q)
    hhg.result[[i]] =  NetworkTest.diffusion.stat(G, X, option = 3, diffusion = TRUE, t.range = 3, n.perm = n.perm, q = diffusion.q)
}

# calculate power using t*
alt.stat = c(); null.stat = c()
pval.mgc = c(); pval.mcorr = c(); pval.hhg = c()
## mgc
for(i in 1:length(mgc.result)){
  if(max(order(-mgc.result[[i]][[1]]))-min(order(-mgc.result[[i]][[1]])) == 2){
    alt.stat[i] = mgc.result[[i]][[1]][mean(order(-mgc.result[[i]][[1]]))]
  }else{
    alt.stat[i] = mgc.result[[i]][[1]][1]
  }

  for(j in 1:nrow(mgc.result[[i]][[2]])){
    if(max(order(-mgc.result[[i]][[2]][j,]))-min(order(-mgc.result[[i]][[2]][j,])) == 2){
      null.stat[j] = mgc.result[[i]][[2]][j,mean(order(-mgc.result[[i]][[2]][j,]))]
    }else{
      null.stat[j] = mgc.result[[i]][[2]][j,1]
    }
  }
  pval.mgc[i] = mean(null.stat >= alt.stat[i])
}

power.mgc = mean(pval.mgc <= 0.05)

## mcorr
for(i in 1:length(mcorr.result)){
  if(max(order(-mcorr.result[[i]][[1]]))-min(order(-mcorr.result[[i]][[1]])) == 2){
    alt.stat[i] = mcorr.result[[i]][[1]][mean(order(-mcorr.result[[i]][[1]]))]
  }else{
    alt.stat[i] = mcorr.result[[i]][[1]][1]
  }

  for(j in 1:nrow(mcorr.result[[i]][[2]])){
    if(max(order(-mcorr.result[[i]][[2]][j,]))-min(order(-mcorr.result[[i]][[2]][j,])) == 2){
      null.stat[j] = mcorr.result[[i]][[2]][j,mean(order(-mcorr.result[[i]][[2]][j,]))]
    }else{
      null.stat[j] = mcorr.result[[i]][[2]][j,1]
    }
  }
  pval.mcorr[i] = mean(null.stat >= alt.stat[i])
}

power.mcorr = mean(pval.mcorr <= 0.05)

## hhg
for(i in 1:length(hhg.result)){
  if(max(order(-hhg.result[[i]][[1]]))-min(order(-hhg.result[[i]][[1]])) == 2){
    alt.stat[i] = hhg.result[[i]][[1]][mean(order(-hhg.result[[i]][[1]]))]
  }else{
    alt.stat[i] = hhg.result[[i]][[1]][1]
  }

  for(j in 1:nrow(hhg.result[[i]][[2]])){
    if(max(order(-hhg.result[[i]][[2]][j,]))-min(order(-hhg.result[[i]][[2]][j,])) == 2){
      null.stat[j] = hhg.result[[i]][[2]][j,mean(order(-hhg.result[[i]][[2]][j,]))]
    }else{
      null.stat[j] = hhg.result[[i]][[2]][j,1]
    }
  }
  pval.hhg[i] = mean(null.stat >= alt.stat[i])
}

power.hhg = mean(pval.hhg <= 0.05)




