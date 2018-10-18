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
#'
library(randomForest)
library(HHG)
library(energy)
library(MGC)
library(copula)
library(dHSIC)
library(MDMR)
library(rerf)
source('GenerateSimulations.R')
test.independence = function(type,n,d,noise,rep){
  if (missing(type)){
    type=6;
  }
  if (missing(n)){
    n=50;
  }
  if (missing(d)){
    d=1;
  }
  if (missing(noise)){
    noise=1;
  }
  if (missing(rep)){
    rep=1000;
  }
  ex=0.5;
  # type=1 to 20, n is suggested to be 100, d=1 with noise=1, or d>1 with noise =0; rep=1000;
  # result=test.independence(type,n,d,noise,rep)
    mgcA=matrix(0,1,rep);
    mgcN=matrix(0,1,rep);
    dcorA=matrix(0,1,rep);
    dcorN=matrix(0,1,rep);
    hhgA=matrix(0,1,rep);
    hhgN=matrix(0,1,rep);
    hsicA=matrix(0,1,rep);
    hsicN=matrix(0,1,rep);
    sRfA=matrix(0,1,rep);
    sRfN=matrix(0,1,rep);
    sRerfA=matrix(0,1,rep);
    sRerfN=matrix(0,1,rep);
    uRfA=matrix(0,1,rep);
    uRfN=matrix(0,1,rep);

    for (i in (1:rep)){
      result=GenerateSimulations(type,n,d,1, noise);
      x=result$x;
      y=result$y;
      dcorA[i] = dcor.ttest(x,y)$estimate;
      hsicA[i] = dhsic(x,y)$dHSIC;
      #Dx = as.matrix(dist(x, diag = TRUE, upper = TRUE))
      Dy = as.matrix(dist(y, diag = TRUE, upper = TRUE))
     # if (ncol(y)>1){
        
      #  forest=RerF(x,y);
      #  sRx=sqrt(1-ComputeSimilarity(x,forest));
      #}
      #else{
       # sRx=sqrt(1-randomForest(x,y,proximity=TRUE)$proximity);
        #forest=RerF(x,y[,1],trees=500);
        #sRx=(1-ComputeSimilarity(x,forest))^ex;
      #}
      Rx=(1-randomForest(x)$proximity)^ex;
      Ry=(1-randomForest(y)$proximity)^ex;
      uRfA[i]=dcor.ttest(Rx,Ry)$estimate;
        
      sRx=(1-randomForest(x,y[,1],proximity=TRUE)$proximity);
      sRfA[i]=dcor.ttest(sRx^ex,Dy)$estimate;
      #sRerfA[i]=dcor.ttest(sRx^0.95,Dy)$estimate;
      #mgcA[i] =  mgc.sample(Dx, Dy)$statMGC;
      #hhgA[i]=hhg.test(Dx,Dy)$sum.chisq;'
      #copulaA[i]=multIndepTest(cbind(x,y),c(ncol(x),ncol(y)),N=11)$global.statistic;
    }

    for (i in (1:rep)){
      result=GenerateSimulations(type,n,d,0, noise);
      x=result$x;
      y=result$y;
      dcorN[i] = dcor.ttest(x,y)$estimate;
      hsicN[i] = dhsic(x,y)$dHSIC;
      #Dx = as.matrix(dist(x, diag = TRUE, upper = TRUE))
      Dy = as.matrix(dist(y, diag = TRUE, upper = TRUE))
      #forest=RerF(x,y[,1],trees=500);
      #sRx=(1-ComputeSimilarity(x,forest))^ex;
      #}
      Rx=(1-randomForest(x)$proximity)^ex;
      Ry=(1-randomForest(y)$proximity)^ex;
      #sRerfN[i]=dcor.ttest(sRx,Dy)$estimate;
      uRfN[i]=dcor.ttest(Rx,Ry)$estimate;
      
      sRx=(1-randomForest(x,y[,1],proximity=TRUE)$proximity);
      sRfN[i]=dcor.ttest(sRx^ex,Dy)$estimate;
      #sRerfN[i]=dcor.ttest(sRx^0.95,Dy)$estimate;
      #mgcN[i] =  mgc.sample(Dx, Dy)$statMGC;
      #hhgN[i]=hhg.test(Dx,Dy,nr.perm=0)$sum.chisq;
      #copulaN[i]=multIndepTest(cbind(x,y),c(ncol(x),ncol(y)),N=11)$global.statistic;
    }

    #testPower.mgc=CalculatePower(mgcA,mgcN);
    testPower.hsic=CalculatePower(hsicA,hsicN);
    testPower.dcor=CalculatePower(dcorA,dcorN);
    #testPower.hhg=CalculatePower(hhgA,hhgN);
    #testPower.copula=CalculatePower(copulaA,copulaN);
    testPower.sRf=CalculatePower(sRfA,sRfN);
    #testPower.sRerf=CalculatePower(sRerfA,sRerfN);
    testPower.uRf=CalculatePower(uRfA,uRfN);
    #testPower.uRf=testPower.sRf;
    #return(list(testPower.mgc=testPower.mgc,testPower.hsic=testPower.hsic,testPower.dcor=testPower.dcor,testPower.hhg=testPower.hhg,testPower.copula=testPower.copula))
    return(list(testPower.dcor=testPower.dcor,testPower.sRf=testPower.sRf,testPower.uRf=testPower.uRf,testPower.hsic=testPower.hsic))
     }

#' This function is an auxliary one that computes testing power
#'
#' @export
#'
CalculatePower = function(A,N,alpha){
  if (missing(alpha)){
    alpha=0.05;
  }
  criticalValue=quantile(N,1-alpha);
  power=mean(A>criticalValue);
}
