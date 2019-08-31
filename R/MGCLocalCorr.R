#' MGC Local Correlations
#'
#' Compute all local correlation coefficients in O(n^2 log n)
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist.X=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist.X=TRUE}.}
#' }
#' @param Y is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{Y is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist.Y=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{Y is a distance matrix. Use flag \code{is.dist.Y=TRUE}.}
#' }
#' @param is.dist.X a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm.X if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params.X a list of trailing arguments to pass to the distance function specified in \code{dist.xfm.X}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return.X the return argument for the specified \code{dist.xfm.X} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm.X[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param is.dist.Y a boolean indicating whether your \code{Y} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm.Y if \code{is.dist == FALSE}, a distance function to transform \code{Y}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{dist.return.Y} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params.Y a list of trailing arguments to pass to the distance function specified in \code{dist.xfm.Y}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return.Y the return argument for the specified \code{dist.xfm.Y} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm.Y(Y)} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm.Y(Y)[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
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
mgc.localcorr <- function(X, Y, is.dist.X=FALSE, dist.xfm.X=mgc.distance, dist.params.X=list(method='euclidean'),
                          dist.return.X=NULL, is.dist.Y=FALSE, dist.xfm.Y=mgc.distance, dist.params.Y=list(method='euclidean'),
                          dist.return.Y=NULL, option='mgc'){
  # validate input is valid and convert to distance matrices, if necessary
  validated <- mgc.validator(X, Y, is.dist.X=is.dist.X, dist.xfm.X=dist.xfm.X, dist.params.X=dist.params.X,
                             dist.return.X=dist.return.X, is.dist.Y=is.dist.Y, dist.xfm.Y=dist.xfm.Y, dist.params.Y=dist.params.Y,
                             dist.return.Y=dist.return.Y)

  DX <- validated$DX; DY <- validated$DY

  return(mgc.localcorr.driver(DX, DY, option=option))
}
#' Driver for MGC Local Correlations
#'
#' @param DX the first distance matrix.
#' @param DY the second distance matrix.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#' @return A list contains the following:
#' \item{\code{corr}}{consists of all local correlations within [-1,1] by double matrix index}
#' \item{\code{varX}}{contains all local variances for X.}
#' \item{\code{varY}}{contains all local variances for X.}
#'
#' @author C. Shen
mgc.localcorr.driver <- function(DX, DY, option='mgc') {

  tmp=mgc.dist.xfm(DX, DY, option)
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

  result=list(corr=corr, varX=varX, varY=varY)
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
#' @keywords internal
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
