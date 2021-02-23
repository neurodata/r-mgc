#' MGC Test
#'
#' The main function that computes the MGC measure between two datasets:
#' It first computes all local correlations, then use the maximal statistic
#' among all local correlations based on thresholding.
#'
#' @references Joshua T. Vogelstein, et al. "Discovering and deciphering relationships across disparate data modalities." eLife (2019).
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
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#'
#' @return A list containing the following:
#' \item{\code{stat}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#' \item{\code{option}}{specifies which global correlation was used}
#'
#' @author C. Shen and Eric Bridgeford
#'
#' @examples
#' library(mgc)
#'
#' n=200; d=2
#' data <- mgc.sims.linear(n, d)
#' mgc.stat.res <- mgc.stat(data$X, data$Y)
#'
#' @export
mgc.stat <- function(X, Y, is.dist.X=FALSE, dist.xfm.X=mgc.distance, dist.params.X=list(method='euclidean'),
                     dist.return.X=NULL, is.dist.Y=FALSE, dist.xfm.Y=mgc.distance, dist.params.Y=list(method='euclidean'),
                     dist.return.Y=NULL, option='mgc') {

  # validate input is valid and convert to distance matrices, if necessary
  validated <- mgc.validator(X, Y, is.dist.X=is.dist.X, dist.xfm.X=dist.xfm.X, dist.params.X=dist.params.X,
                             dist.return.X=dist.return.X, is.dist.Y=is.dist.Y, dist.xfm.Y=dist.xfm.Y, dist.params.Y=dist.params.Y,
                             dist.return.Y=dist.return.Y)

  DX <- validated$DX; DY <- validated$DY

  return(mgc.stat.driver(DX, DY, option=option))
}

#' MGC Sample Statistic Internal Driver
#' @param DX the first distance matrix.
#' @param DY the second distance matrix.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#' @keywords internal
mgc.stat.driver <- function(DX, DY, option='mgc') {
  # compute local correlation map
  localCorr <- mgc.localcorr.driver(DX, DY, option)$corr # compute all localCorr
  m <- nrow(localCorr)
  n <- ncol(localCorr)

  if (m==1 || n==1){
    stat <- localCorr[m, n]
    optimalScale=m*n
  } else {
    sz <- nrow(DX) - 1 # sample size minus 1
    R <- Thresholding(localCorr, m, n, sz) # find a connected region of significant local correlations
    res <- Smoothing(localCorr,m, n, R) # find the maximal within the significant region
    optimalScale <- res$optimalScale; stat <- res$stat
  }

  return(list(stat=stat, localCorr=localCorr, optimalScale=optimalScale, option=option))
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
#' @author Eric Bridgeford and C. Shen
#' @keywords internal
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
#' \item{\code{stat}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{optimalScale}}{the estimated optimal scale as a list.}
#'
#' @author C. Shen
#' @keywords internal
Smoothing <- function(localCorr,m,n,R){
  stat=localCorr[m,n] # default sample mgc to local corr at maximal scale
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
        if (tmp>=stat){
          stat=tmp
          optimalScale=list(x=k, y=l) # take the scale of maximal stat and change to single index
        }
      #}
    }
  }
  result=list(stat=stat,optimalScale=optimalScale)
  return(result)
}

#' Connected Components Labelling -- Unique Patch Labelling
#'
#' \code{ConnCompLabel} is a 1 pass implementation of connected components
#' labelling. Here it is applied to identify disjunt patches within a
#' distribution. \cr \cr The raster matrix can be a raster of class 'asc'
#' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package).
#'
#'
#' @param mat is a binary matrix of data with 0 representing background and 1
#' representing environment of interest. NA values are acceptable. The matrix
#' can be a raster of class 'asc' (this & adehabitat package), 'RasterLayer'
#' (raster package) or 'SpatialGridDataFrame' (sp package)
#' @return A matrix of the same dim and class of \code{mat} in which unique
#' components (individual patches) are numbered 1:n with 0 remaining background
#' value.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references Chang, F., C.-J. Chen, and C.-J. Lu. 2004. A linear-time
#' component-labeling algorithm using contour tracing technique. Comput. Vis.
#' Image Underst. 93:206-220.
#' @examples
#'
#'
#' #define a simple binary matrix
#' tmat = { matrix(c( 0,0,0,1,0,0,1,1,0,1,
#'                    0,0,1,0,1,0,0,0,0,0,
#'                    0,1,NA,1,0,1,0,0,0,1,
#'                    1,0,1,1,1,0,1,0,0,1,
#'                    0,1,0,1,0,1,0,0,0,1,
#'                    0,0,1,0,1,0,0,1,1,0,
#'                    1,0,0,1,0,0,1,0,0,1,
#'                    0,1,0,0,0,1,0,0,0,1,
#'                    0,0,1,1,1,0,0,0,0,1,
#'                    1,1,1,0,0,0,0,0,0,1),nr=10,byrow=TRUE) }
#'
#' #do the connected component labelling
#' ccl.mat = ConnCompLabel(tmat)
#' ccl.mat
#' image(t(ccl.mat[10:1,]),col=c('grey',rainbow(length(unique(ccl.mat))-1)))
#'
#'
#' @useDynLib mgc
#' @importFrom raster setValues
#' @export
ConnCompLabel <- function(mat)	{
  attrib = attributes(mat)
  #check to ensure matrix
  mat = try(as.matrix(mat))
  if (!is.matrix(mat)) stop('objects must be a matrix')
  #run the connected component labelling
  out = .Call('ccl',mat,PACKAGE='mgc')
  #reset the attributes of the input
  if (any(class(attrib) %in% 'RasterLayer')) {
    attrib = setValues(attrib, as.vector(t(t(unclass(out))[dim(out)[2]:1,]))); return(attrib)
  } else if (any(class(attrib) == 'SpatialGridDataFrame')) {
    attrib@data[1] = as.vector(unclass(out)[,dim(out)[2]:1]); return(attrib)
  } else {
    attributes(out) = attrib; return(out)
  }
}
