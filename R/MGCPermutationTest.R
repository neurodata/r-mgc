#' The main function that tests independent between two data sets by MGC and permutation test.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a [n x n] distance matrix}{X is a square matrix with zeros on diagonal}
#'    \item{a [n x d] data matrix}{Otherwise}
#' }
#' @param Y is interpreted as:
#' \describe{
#'    \item{a [n x n] distance matrix}{Y is a square matrix with zeros on diagonal}
#'    \item{a [n x d] data matrix}{Otherwise}
#' }
#' @param rep=1000 specifies the number of replicates to use for the permutation test;
#' @param option='mgc' is a string that specifies which global correlation to build up-on.
#' \describe{
#'    \item{'mgc'}{use the MGC global correlation.}
#'    \item{'dcor'}{use the dcor global correlation.}
#'    \item{'mantel'}{use the mantel global correlation.}
#'    \item{'rank'}{use the rank global correlation.}
#' }
#' @return P-value and test statistic of MGC;
#' @return All local p-values by double matrix index, all local correlations by double matrix index, and the estimated optimal scales as matrix single indices.
#'
#' Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
#' unless corrected for multiple testing problem.
#'
#' @author C. Shen
#' @export
#'
mgc.test <-function(X, Y, rep=1000, option='mgc'){
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(as.matrix(X)) != ncol(as.matrix(X)) | sum(diag(X)^2) > 0){
    X = as.matrix(dist(X, method='euclidean'))
    # print('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
  }
  if (nrow(as.matrix(Y))!=ncol(as.matrix(Y))|sum(diag(Y)^2)>0){
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
  for (r in (1:rep)){
    # Use random permutations on the second data set
    per = sample(np)
    YN = Y[per, per]
    tmp = mgc.sample(X, YN, option)
    pMGC = pMGC + (tmp$statMGC >= result$statMGC)*1/rep
    pLocalCorr = pLocalCorr + (tmp$localCorr >= result$localCorr)*1/rep
  }

  # Estimate the optimal scales via both the statistics and p-values if possible
  # optimalScale=FindLargestRectangles((pLocalCorr<=pMGC)&(result$localCorr>=result$statMGC))$M;
  # optimalScale=which(optimalScale==1);
  # if (length(optimalScale)==0){
  #  optimalScale=FindLargestRectangles((pLocalCorr<=pMGC))$M; # If empty, estimate the optimal scales via the p-value only
  #  optimalScale=which(optimalScale==1);
  #  if (length(optimalScale)==0){
  #    optimalScale=m*n; # default
  #  }
  #}

  result=list(pMGC=pMGC,statMGC=result$statMGC,pLocalCorr=pLocalCorr,localCorr=result$localCorr,optimalScale=result$optimalScale)
  return(result)
}
