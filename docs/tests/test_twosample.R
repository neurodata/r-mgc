#' This function compute the testing power for two-sample testing, based on the two
#' dimension simulations with rotation matrix.
#'
#' @param type is a number between 1 to 20,
#' @param n is the sample size,
#' @param deg specifies the rotation angle, say 0 to 2\pi
#' @param noise specifies the amount of noise
#' @param rep is the number of replicates for estimating power
#' @return result: list of testing powers, for mgc/hsic/dcor/hhg
#' @export
#'

library(HHG)
library(energy)
library(dHSIC)
source('GenerateSimulations.R')

test.twosample = function(type,n,deg,noise,rep){
  if (missing(type)){
    type=6;
  }
  if (missing(n)){
    n=50;
  }
  if (missing(deg)){
    deg=pi/8;
  }
  if (missing(noise)){
    noise=0;
  }
  if (missing(rep)){
    rep=100;
  }
  alpha=0.05;
  d=1;
  rotate=matrix(c(cos(deg),sin(deg),-sin(deg),cos(deg)),2,2);
  pmgc=matrix(0,1,rep);
  pdcor=matrix(0,1,rep);
  phhg=matrix(0,1,rep);
  phsic=matrix(0,1,rep);

  rep2=100; # replicates for permutation test
  for (i in (1:rep)){
    result=GenerateSimulations(type,n,d,1, noise);
    x1=cbind(result$x,result$y);
    result=GenerateSimulations(type,n,d,1, noise);
    x2=cbind(result$x,result$y);
    x=rbind(x1,x2%*%rotate);
    y=rbind(matrix(0,n,1),matrix(1,n,1));
    Dx = as.matrix(dist(x, diag = TRUE, upper = TRUE))
    Dy = as.matrix(dist(y, diag = TRUE, upper = TRUE))

    phsic[i] = dhsic.test(x,y,method = "permutation",B=rep2)$p.value;
    pdcor[i] = dcov.test(x,y,R=rep2)$p.value;
    pmgc[i] =  mgc.test(Dx, Dy,rep2)$pMGC;
    phhg[i]=hhg.test(Dx,Dy,nr.perm=rep2)$perm.pval.hhg.sc;
  }

  testPower.mgc=mean(pmgc<alpha)
  testPower.hsic=mean(phsic<alpha)
  testPower.dcor=mean(pdcor<alpha)
  testPower.hhg=mean(phhg<alpha)
  return(list(testPower.mgc=testPower.mgc,testPower.hsic=testPower.hsic,testPower.dcor=testPower.dcor,testPower.hhg=testPower.hhg))
}
