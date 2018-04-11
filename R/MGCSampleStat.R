#' MGC Sample
#'
#' The main function that computes the MGC measure between two datasets:
#' It first computes all local correlations,
#' then use the maximal statistic among all local correlations based on thresholding.
#'
#' @param A is interpreted as:
#' \describe{
#'    \item{a \code{[n x n]} distance matrix}{A is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{A is a data matrix with \code{n} samples in \code{d} dimensions.}
#' }
#' @param B is interpreted as:
#' \describe{
#'    \item{a \code{[n x n]} distance matrix}{B is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{B is a data matrix with \code{n} samples in \code{d} dimensions.}
#' }
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#' @return A list containing the following:
#' \item{\code{statMGC}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{localCorr}}{consists of all local correlations by double matrix index}
#' \item{\code{optimalScale}}{the estimated optimal scale in matrix single index.}
#' @author C. Shen
#'
#' @examples
#' library(mgc)
#'
#' n=200; d=2
#' data <- mgc.sims.linear(n, d)
#' lcor <- mgc.sample(data$X, data$Y)
#'
#' @export
mgc.sample <- function(A, B, option='mgc'){
  localCorr=mgc.localcorr(A,B,option)$corr # compute all localCorr
  m=nrow(localCorr)
  n=ncol(localCorr)
  if (m==1||n==1){
    statMGC=localCorr[m,n]
    optimalScale=m*n
  } else {
    sz=nrow(A)-1 # sample size minus 1
    R=Thresholding(localCorr,m,n,sz) # find a connected region of significant local correlations
    res=Smoothing(localCorr,m,n,R) # find the maximal within the significant region
  }

  result=list(statMGC=res$statMGC,localCorr=localCorr,optimalScale=res$optimalScale)
  return(result)
}

#' An auxiliary function that finds a region of significance in the local correlation map by thresholding.
#'
#' @importFrom stats qbeta
#'
#' @param localCorr is all local correlations
#' @param m is the number of rows of localCorr
#' @param n is the number of columns of localCorr
#' @param sz is the sample size of original data (which may not equal m or n in case of repeating data).
#'
#' @return R is a binary matrix of size m and n, with 1's indicating the significant region.
#' @author C. Shen
#' @importFrom SDMTools ConnCompLabel
#'
Thresholding <- function(localCorr,m,n,sz){
  # A threshold is estimated based on normal distribution approximation from Szekely2013
  prt=1-0.02/sz # percentile to consider as significant
  # thres=sqrt((sz*(sz-3)/2)-1) # normal approximation, which is equivalent to beta approximation for n larger than 10
  # thres=qnorm(prt)/thres
  thres=sz*(sz-3)/4-1/2 # beta approximation
  thres=(qbeta(prt,thres,thres))*2-1

  opt=0 # set opt=1 and add the negative local correlation as a non-parametric and data-adaptive threshold
  if (opt==1){
    thres1=localCorr
    thres1=thres1[thres1<0] # all negative correlations
    thres1=5*sqrt(sum(thres1^2)/length(thres1))  # the standard deviation of negative correlations
    # Use the maximal of paratemetric and non-parametric thresholds
    if (is.na(thres1)==FALSE && thres1>thres){
      thres=thres1
    }
  }
  thres=max(thres,localCorr[m,n]) # take the maximal of threshold and local correlation at the maximal scale

  # Find the largest connected component of significant correlations
  R=(localCorr>thres)
  if (sum(R)>0){
    R=ConnCompLabel(R==1)
    tmp=tabulate(R)
    tmp=which.max(tmp)
    R=(R==tmp)
  } else {
    R=0
  }
  return(as.matrix(R))
}

#' An auxiliary function that finds the smoothed maximal within the significant region R:
#' If area of R is too small, return the last local corr otherwise take the maximum within R.
#'
#' @param localCorr is all local correlations
#' @param m is the number of rows of localCorr
#' @param n is the number of columns of localCorr
#' @param R is a binary matrix of size m by n indicating the significant region.
#'
#' @return A list contains the following:
#' \item{\code{statMGC}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{optimalScale}}{the estimated optimal scale as a list.}
#'
#' @author C. Shen
Smoothing <- function(localCorr,m,n,R){
  statMGC=localCorr[m,n] # default sample mgc to local corr at maximal scale
  optimalScale=list(x=m, y=n) # default the optimal scale to maximal scale
  if (norm(R,"F")!=0){
    # tau=1 # number of adjacent scales to smooth with
    if (sum(R)>=2*(min(m,n))){ # proceed only when the region area is sufficiently large
      tmp=max(localCorr[R])
      ind=which((localCorr>=tmp)&(R==1)) # find all scales within R that maximize the local correlation
      k = ((ind-1) %% m) + 1
      l = floor((ind-1) / m) + 1

      #ln=ceiling(tau) # number of adjacent rows to check
      #km=ceiling(tau) # number of adjacent columns to check
      #for (i in (1:length(k))){
       # ki=k[i]
      #  li=l[i]

        # index of adjacent rows and columns
       # left=max(2,li-ln)
      #  right=min(n,li+ln)
       # upper=max(2,ki-km)
      #  down=min(m,ki+km)

       # tmp1=min(localCorr[upper:down,li]) # take minimal correlation at given row and along adjacent columns
      #  tmp2=min(localCorr[ki,left:right]) # take minimal correlation at given column and along adjacent rows
       # tmp=max(tmp1,tmp2) # take the min for sample mgc
        if (tmp>=statMGC){
          statMGC=tmp
          optimalScale=list(x=k, y=l) # take the scale of maximal stat and change to single index
        }
      #}
    }
  }
  result=list(statMGC=statMGC,optimalScale=optimalScale)
  return(result)
}
