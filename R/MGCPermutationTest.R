#' The main function that tests independent between two data sets by MGC and permutation test.
#'
#' @param X is a distance matrix or a n*d data matrix; if it is not a square matrix with zeros on diagonal, it is treated as n*d data; 
#' @param Y is a second distance matrix or a n*d data matrix, with the same distance matrix check as X;
#' @param rep specifies the number of replicates to use for the permutation test;
#' @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank'.
#'
#' @return A list contains the following output:
#' @return P-value and test statistic of MGC;
#' @return All local p-values by double matrix index, all local correlations by double matrix index, and the estimated optimal scales as matrix single indices.
#'
#' Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
#' unless corrected for multiple testing problem.
#'
#' @export
#' 
MGCPermutationTest <-function(X,Y,rep,option){
  if (missing(rep)){
    rep=1000; # use 1000 random permutations by default
  }
  if (missing(option)){
    option='mgc'; # use mgc by default
  }
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(as.matrix(X))!=ncol(as.matrix(X))|sum(diag(X)^2)>0){
    X=as.matrix(dist(X,method='euclidean'));
    # print('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
  }
  if (nrow(as.matrix(Y))!=ncol(as.matrix(Y))|sum(diag(Y)^2)>0){
    Y=as.matrix(dist(Y,method='euclidean'));
    # print('The second data is not a Euclidean distance matrix; transformed to distance matrix instead.')
  }
  np=nrow(Y);

  # Compute sample MGC and all local correlations
  result=MGCSampleStat(X,Y,option);
  m=nrow(result$localCorr);
  n=ncol(result$localCorr);
  pLocalCorr=matrix(0,m,n);
  pMGC=0;

  # Compute sample MGC and all local correlations for each permuted data
  for (r in (1:rep)){
    # Use random permutations on the second data set
    per=sample(np);
    YN=Y[per,per];
    tmp=MGCSampleStat(X,YN,option);
    pMGC=pMGC+(tmp$statMGC>=result$statMGC)*1/rep;
    pLocalCorr=pLocalCorr+(tmp$localCorr>=result$localCorr)*1/rep;
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

  result=list(pMGC=pMGC,statMGC=result$statMGC,pLocalCorr=pLocalCorr,localCorr=result$localCorr,optimalScale=result$optimalScale);
  return(result);
}
