#' MGC Local Correlations
#'
#' Compute all local correlation coefficients in O(n^2 log n)
#'
#' @param X is interpreted as:
#' \itemize{
#'    \item{a \code{[n x n]} distance matrix}{X is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions.}
#'  }
#' @param Y is interpreted as:
#' \describe{
#'    \item{a \code{[n x n]} distance matrix}{Y is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{Y is a data matrix with \code{n} samples in \code{d} dimensions.}
#' }
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{'mgc'}{use the MGC global correlation.}
#'    \item{'dcor'}{use the dcor global correlation.}
#'    \item{'mantel'}{use the mantel global correlation.}
#'    \item{'rank'}{use the rank global correlation.}
#' }
#' @return A list contains the following:
#' \item{\code{corr}}{consists of all local correlations within [-1,1] by double matrix index}
#' \item{\code{varX}}{contains all local variances for X.}
#' \item{\code{varY}}{contains all local variances for X.}
#'
#' @author C. Shen
#'
#' @examples
#' library(mgc)
#'
#' n=200; d=2
#' data <- mgc.sims.linear(n, d)
#' lcor <- mgc.localcorr(data$X, data$Y)
#'
#' @export
#'
mgc.localcorr <- function(X,Y,option='mgc'){
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(as.matrix(X))!=ncol(as.matrix(X))|sum(diag(X)^2)>0){
    X=as.matrix(dist(X,method='euclidean'))
    # print('The first data is not a Euclidean distance matrix transformed to distance matrix instead.')
  }
  if (nrow(as.matrix(Y))!=ncol(as.matrix(Y))|sum(diag(Y)^2)>0){
    Y=as.matrix(dist(Y,method='euclidean'))
    # print('The second data is not a Euclidean distance matrix transformed to distance matrix instead.')
  }

  tmp=mgc.distTransform(X,Y,option)
  corr=LocalCov(tmp$A,t(tmp$B),tmp$RX,t(tmp$RY)) # compute all local covariances
  varX=LocalCov(tmp$A,t(tmp$A),tmp$RX,t(tmp$RX)) # compute local variances for first data
  varY=LocalCov(tmp$B,t(tmp$B),tmp$RY,t(tmp$RY)) # compute local variances for second data
  varX=diag(varX)
  varY=diag(varY)

  options(warn=-1)
  corr=corr/Re(sqrt(varX%*%t(varY)))
  corr[corr>1]=1 # avoid computational issue that may cause a few local corr to be negligably larger than 1
  options(warn=0)

  # set any local correlation to 0 if any corresponding local variance is no larger than 0
  for (k in (1:length(varX))){
    if (varX[k]<=0){
      corr[k,]=0
    }
  }
  for (l in (1:length(varY))){
    if (varY[l]<=0){
      corr[,l]=0
    }
  }

  result=list(corr=corr,varX=varX,varY=varY)
  return(result)
}

#' An auxiliary function that computes all local correlations simultaneously in O(n^2).
#'
#' @param A is a properly transformed distance matrix
#' @param B is the second distance matrix properly transformed
#' @param RX is the column-ranking matrix of A
#' @param RY is the column-ranking matrix of B.
#'
#' @return covXY is all local covariances computed iteratively.
#'
LocalCov <- function(A,B,RX,RY){
  n=nrow(A); nX=max(RX); nY=max(RY)
  covXY=matrix(0,nX,nY)#varX=rep(0,nX)varY=rep(0,nY)
  EX=rep(0,nX); EY=rep(0,nY)

  # summing up the entriwise product of A and B based on the ranks, which
  # yields the local family of covariance and variances
  for (j in (1:n)){
    for (i in (1:n)){
      a=A[i,j]
      b=B[i,j]
      k=RX[i,j]
      l=RY[i,j]
      covXY[k,l]=covXY[k, l]+a*b
      # varX[k]=varX[k]+a^2
      # varY[l]=varY[l]+b^2
      EX[k]=EX[k]+a
      EY[l]=EY[l]+b
    }
  }

  for (k in (1:(nX-1))){
    covXY[k+1,1]=covXY[k,1]+covXY[k+1,1]
    # varX[k+1]=varX[k]+varX[k+1]
    EX[k+1]=EX[k]+EX[k+1]
  }
  for (l in (1:(nY-1))){
    covXY[1,l+1]=covXY[1,l]+covXY[1,l+1]
    # varY[l+1]=varY[l]+varY[l+1]
    EY[l+1]=EY[l]+EY[l+1]
  }
  for (l in (1:(nY-1))){
    for (k in (1:(nX-1))){
      covXY[k+1,l+1]=covXY[k+1,l]+covXY[k,l+1]+covXY[k+1,l+1]-covXY[k,l]
    }
  }

  # normalize the covariance by the variances yields the local family of correlation
  covXY=(covXY-EX%*%t(EY)/n/(n))
  # varX=varX-EX^2/n^2
  # varY=varY-EY^2/n^2
  # covXY[1,1:nY]=0 # local cov without any neighbor is meaningless and set to 0 instead
  # covXY[1:nX,1]=0

  return(covXY)
}

#LocalWeights <- function(A,B,RX,RY,ind){
# An auxiliary function that computes the contributions of each distance entries to
# the local distance correlation at a given scale.
# nX=max(RX)nY=max(RY)
#if (ind>nX*nY || ind<1){
# ind=nX*nY # default to global scale when the specified index is out of range
#}
#k = ((ind-1) %% nX) + 1
#l = floor((ind-1) / nX) + 1
#RX=(RX>k)
#RY=(RY>l)
#A[RX]=0
#B[RY]=0
#weight=(A-mean(A))*(B-mean(B))
#return(weight)
#}
