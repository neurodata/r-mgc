#setwd("/users/ylee/multi")

#source("/users/ylee/newMGC/MGCDistTransform.R")
#source("/users/ylee/newMGC/FindLargestRectangles.R")
#source("/users/ylee/newMGC/MGCLocalCorr.R")
#source("/users/ylee/newMGC/MGCSampleStat.R")
#source("/users/ylee/newMGC/MGCPermutationTest.R")

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
library(nparcomp)
library(Hotelling)
library(randomForest)
library(MGC)

n = 50 # sample size
nr = 100 # number of independent iterations
epsilon = seq(0, pi/2,  pi/10)
power = matrix(0, length(epsilon),5)

for(i in 1:nr){
  print(i)
  set.seed(i)
  mgc = c(); mcor = c(); hhg = c(); hotel = c(); energy = c()

  for(k in 1:length(epsilon)){

    affine = matrix(c(cos(epsilon[k]), sin(epsilon[k]), -sin(epsilon[k]), cos(epsilon[k])), 2, 2)

    X1 = runif(n, -1, 1)
    X2 = X1 + rnorm(n, 0, 1)
    X = cbind(X1, X2)

    tmp.Y1 = runif(n, -1, 1)
    tmp.Y2 = tmp.Y1 + rnorm(n, 0, 1)
    tmp.Y = cbind(tmp.Y1, tmp.Y2)
    Y = tmp.Y %*% affine

    Dz = as.matrix(dist(rbind(X,Y), diag = TRUE, upper = TRUE))
    L = c(rep(0, nrow(X)), rep(1, nrow(Y)))
    Dl = as.matrix(dist(L, diag = TRUE, upper = TRUE))


    ## distance-based test
    tmp =  MGCPermutationTest(Dz, Dl,100, option = 'mgc')
    mgc[k] = tmp[[1]]
    mcor[k] = tmp[[3]][dim(tmp[[3]])[1],dim(tmp[[3]])[2]]
    hhg[k] = hhg.test(Dz, Dl, nr.perm = 100)$perm.pval.hhg.sl
    hotel[k] = hotelling.test(X, Y)$pval

    ## distance-based test
    energy[k] = disco(rbind(X,Y), factors = as.factor(L), method = "discoB",
    index=1.0, R=99)$p.value
  }
   pvals = cbind(mgc, mcor, hhg, hotel, energy)
   power=power+ as.matrix(pvals <= 0.05)/nr;
}

