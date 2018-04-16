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
    rep=1000;
  }
  n=seq(5,100,5);

  sz=20;
  dcor=matrix(0,1,20);
  copula=matrix(0,1,20);
  hsic=matrix(0,1,20);

  #for (i in (1:20)){
    for (j in (1:sz)){
    power=test.independence(type,n[j],1,1,rep);
    dcor[j]=power$testPower.dcor;
    hsic[j]=power$testPower.hsic;
    copula[j]=power$testPower.copula;
  #}
  }
  #return(list(testPower.mgc=testPower.mgc,testPower.hsic=testPower.hsic,testPower.dcor=testPower.dcor,testPower.hhg=testPower.hhg,testPower.copula=testPower.copula))
  return(list(testPower.dcor=dcor,testPower.copula=copula,testPower.hsic=hsic,n=n))
}
