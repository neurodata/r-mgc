#' Transform the distance matrices, with column-wise ranking if needed.
#'
#' @param X is a distance matrix;
#' @param Y is a second distance matrix;
#' @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank';
#' @param optionRk is a string that specifies whether ranking within column is computed or not.
#'
#' @return A list contains the following output:
#' @return A and B are the centered distance matrices;
#' @return RX and RY are the column rank matrices of X and Y respectively.
#'
#' @export
#'
MGCDistTransform <- function(X,Y,option,optionRk){
  if (missing(option)){
    option='mgc'; # default to mgc transform
  }
  if (missing(optionRk) || option=='rank'){
    optionRk=1; # do ranking or not, 0 to no ranking
  }

  # Depending on the choice of the global correlation, properly transform each distance matrix
  tA=DistCentering(X,option,optionRk);
  tB=DistCentering(Y,option,optionRk);

  result=list(A=tA$A,RX=tA$RX,B=tB$A,RY=tB$RX);
  return(result);
}

#' An auxiliary function that properly transforms the distance matrix X
#'
#' @param X is a symmetric distance matrix;
#' @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank';
#' @param optionRk is a string that specifies whether ranking within column is computed or not.
#'
#' @return A list contains the following output:
#' @return A is the centered distance matrices;
#' @return RX is the column rank matrices of X.
#'
DistCentering<-function(X,option,optionRk){
  n=nrow(X);
  if (optionRk!=0){
    RX=DistRanks(X); # the column ranks for X
  } else {
    RX=matrix(0,n,n);
  }

  if (option=='rank'){
    X=RX;
  }
  # Default mgc transform
  EX=t(matrix(rep(colMeans(X)*n/(n-1),n), ncol = n));

  if (option=='dcor'){ # unbiased dcor transform
    EX=t(matrix(rep(colMeans(X)*n/(n-2),n), ncol = n))+matrix(rep(rowMeans(X)*n/(n-2),n), ncol = n)-sum(X)/(n-1)/(n-2);
    EX=EX+X/n;
  }
  if (option=='mantel'){ # mantel transform
    EX=sum(X)/n/(n-1);
  }
  # if (option=='dcorDouble'){ # original double centering of dcor
  # EX=t(matrix(rep(colMeans(X),n), ncol = n))+matrix(rep(rowMeans(X),n), ncol = n)-mean(X);
  # }
  A=X-EX;

  # The diagonal entries are always excluded
  for (j in (1:n)){
    A[j,j]=0;
  }

  result=list(A=A,RX=RX);
  return(result);
}

#' An auxiliary function that sorts the entries within each column by ascending order:
#' For ties, the minimum ranking is used,
#' e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.
#'
#' @param dis is a symmetric distance matrix.
#'
#' @return disRank is the column rank matrices of X.
#'
DistRanks <- function(dis) {
  n=nrow(dis);
  disRank=matrix(0,n,n);
  for (i in (1:n)){
    v=dis[,i];
    tmp=rank(v,ties.method="min");
    tu=unique(tmp);
    if (length(tu)!=max(tmp)){
      tu=sort(tu);
      for (j in (1:length(tu))){
        kk=which(tmp==tu[j]);
        tmp[kk]=j;
      }
    }
    disRank[,i]=tmp;
  }
  return(disRank);
}
