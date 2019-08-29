#' This function compute the testing power for independence testing, based on the two
#' dimension simulations with rotation matrix.
#'
#' @param type is a number between 1 to 20,
#' @param n is the sample size,
#' @param d is the dimensionality that is an integer no smaller than 1.
#' @param noise specifies the amount of noise
#' @param rep is the number of replicates for estimating power
#' @return result: list of testing powers, for mgc/hsic/dcor/hhg
#' @export
#'
#'source('run_1d_sims.R')

source('test_independence.R')
run1DSims = function(type,rep){
  if (missing(rep)){
    rep=100;
  }
  sz=10; # reduce sz for faster running
  n=seq(100/sz,100,100/sz);

  dcor=matrix(0,1,sz);
  uRf=matrix(0,1,sz);
  sRf=matrix(0,1,sz);
  #sRerf=matrix(0,1,sz);
  hsic=matrix(0,1,sz);

  #for (i in (1:20)){
    for (j in (1:sz)){
    power=test.independence(type,n[j],1,1,rep);
    dcor[j]=power$testPower.dcor;
    hsic[j]=power$testPower.hsic;
    uRf[j]=power$testPower.uRf;
    sRf[j]=power$testPower.sRf;
    #sRerf[j]=power$testPower.sRerf;
  #}
    }
  
  avePower=cbind(mean(sRf),mean(uRf),mean(dcor),mean(hsic));
  #return(list(testPower.mgc=testPower.mgc,testPower.hsic=testPower.hsic,testPower.dcor=testPower.dcor,testPower.hhg=testPower.hhg,testPower.copula=testPower.copula))
  return(list(average=avePower,testPower.dcor=dcor,testPower.sRf=sRf,testPower.uRf=uRf,testPower.hsic=hsic,n=n))
}
