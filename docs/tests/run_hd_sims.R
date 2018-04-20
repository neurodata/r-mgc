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
#'source('run_hd_sims.R')

source('test_independence.R')
runHDSims = function(type,rep){
  if (missing(rep)){
    rep=1000;
  }
  n=100;
  sz=20;
  dcor=matrix(0,1,sz+1);
  copula=matrix(0,1,sz+1);
  hsic=matrix(0,1,sz+1);

  #for (i in (1:20)){
    dimUpper=fixDim(type);
    sz=20;
    if (dimUpper<20){
      sz=dimUpper;
    }
    if (dimUpper/sz==1){
      dim=seq(1,dimUpper,1);
    }else{
      dim=seq(ceiling(dimUpper/sz),dimUpper,ceiling(dimUpper/sz));
      dim=c(1,dim);
      sz=sz+1;
    }
    for (j in (1:sz)){
      if (dim[j]!=0){
      power=test.independence(type,n,dim[j],0,rep);
      dcor[j]=power$testPower.dcor;
      hsic[j]=power$testPower.hsic;
      copula[j]=power$testPower.copula;
      }
    }
  #}
  #return(list(testPower.mgc=testPower.mgc,testPower.hsic=testPower.hsic,testPower.dcor=testPower.dcor,testPower.hhg=testPower.hhg,testPower.copula=testPower.copula))
  return(list(testPower.dcor=dcor,testPower.copula=copula,testPower.hsic=hsic,dim=dim));
}

fixDim=function(type){
  dimUpper=20;
  if ((type==1) || (type==2) || (type==3)){
    dimUpper=1000;
  }
  if ((type==4) || (type==12) || (type==13) || (type==19)){
    dimUpper=10;
  }
  if ((type==14) || (type==18)) {
    dimUpper=40;
  }
  if ((type==9) || (type==10) || (type==20)){
    dimUpper=100;
  }
  return(dimUpper);
}
