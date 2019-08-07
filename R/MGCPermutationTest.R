#' MGC Permutation Test
#'
#' The main function that tests independent between two data sets by MGC and permutation test.
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
#' as the \code{$D} return argument. See \link[mgc]{discr.distance} for details.
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
#' as the \code{dist.return.Y} return argument. See \link[mgc]{discr.distance} for details.
#' @param dist.params.Y a list of trailing arguments to pass to the distance function specified in \code{dist.xfm.Y}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return.Y the return argument for the specified \code{dist.xfm.Y} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm.Y(Y)} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm.Y(Y)[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param nperm specifies the number of replicates to use for the permutation test. Defaults to \code{1000}.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#'
#' @return A list containing the following:
#' \item{\code{p.value}}{P-value of MGC}
#' \item{\code{stat}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{pLocalCorr}}{P-value of the local correlations by double matrix index}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#'
#' Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
#' unless corrected for multiple testing problem.
#'
#' @author C. Shen
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("mgc", package = "mgc")}
#'
#' @examples
#'
#' library(mgc)
#'
#' n = 100; d = 2
#' data <- mgc.sims.linear(n, d)
#' result <- mgc.test(data$X, data$Y, nperm=10)
#'
#' @export
mgc.test <-function(X, Y, nperm=1000, option='mgc'){
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(as.matrix(X)) != ncol(as.matrix(X)) | sum(diag(X)^2) > 0){
    X = as.matrix(dist(X, method='euclidean'))
    # print('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
  }
  if (nrow(as.matrix(Y)) != ncol(as.matrix(Y)) | sum(diag(Y)^2)>0){
    Y=as.matrix(dist(Y, method='euclidean'))
    # print('The second data is not a Euclidean distance matrix; transformed to distance matrix instead.')
  }
  np = nrow(Y)

  # Compute sample MGC and all local correlations
  result = mgc.sample(X, Y, option)
  m = nrow(result$localCorr)
  n = ncol(result$localCorr)
  pLocalCorr = matrix(0,m,n)
  pMGC = 0

  # Compute sample MGC and all local correlations for each permuted data
  for (r in (1:nperm)){
    # Use random permutations on the second data set
    per = sample(np)
    YN = Y[per, per]
    tmp = mgc.sample(X, YN, option)
    pMGC = pMGC + (tmp$statMGC >= result$statMGC)*1/nperm
    pLocalCorr = pLocalCorr + (tmp$localCorr >= result$localCorr)*1/nperm
  }

  result=list(p.value=pMGC,stat=result$statMGC,pLocalCorr=pLocalCorr,localCorr=result$localCorr,optimalScale=result$optimalScale)
  return(result)
}
