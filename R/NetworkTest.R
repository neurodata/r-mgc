#' This function prints out the distance-based independence test (MGC, mCorr, HHG) statistics.
#' The statistics depend on the choice of diffusion map embeddings. 
#'
#' @param G is a igraph object having n nodes;
#' @param X is a vector of length [n] or [nxn] matrix;
#' @param option specifies the types of distance-based independence test statistic.
#' \describe{
#'    \item{'1'}{use the MGC.}
#'    \item{'2'}{use the mCorr.}
#'    \item{'3'}{use the HHG.}
#' }
#' @importFrom igraph get.adjacency
#' @param diffusion is interpreted as:
#' \describe{
#'    \item{'TRUE'}{use the diffusion map embeddings.}
#'    \item{'FALSE'}{use the columns of adjacency matrix as the ad-hoc graph embeddings.}
#' }
#' @param t.range is a range of Markov iteration times applied for diffusion map embeddings;
#' @param n.perm is the number of permutation samples.
#' @return p-value of the test if \code{diffusion} is \code{"FALSE"};
#' @return a collection of statistics under the alternative and \code{n.perm} number of nulls across \code{t.range} if \code{diffusion} is \code{"TRUE"}.
#' @author Youjin Lee
#' @export
#'

NetworkTest.diffusion.stat = function(G, X, option, diffusion, t.range, n.perm){
  
  A = as.matrix(get.adjacency(G))
  D = diag(pmin( (rowSums(A))^(-1/2) , 1))
  P = D %*% A %*% D # a nomalized graph Laplacian
  
  if(diffusion == FALSE){ 
    Dx = as.matrix(dist((A)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      pvalues = mgc.test(Dx, Dy, rep = n.perm, option = 'mgc')$pMGC
      return(pvalues)
    }else if(option == 2){
      pvalues = dcov.test(Dx, Dy, index = 1.0, R = n.perm)$p.value
      return(pvalues)
    }else if(option == 3){
      pvalues = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
      return(pvalues)
    }
  }
  
  mgc.stat = c(); dcor.stat = c(); hhg.stat = c()
  
  for(s in 1:length(t.range)){

    # as a dimensional choice, use the second elbow of absolute eivengalues from a diffusion map at t = 1.
    diffusion.q  =  min(max(getElbows(abs(print.lambda(P, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A)-1)
    U  =  dmap.q(P, t.range[s], diffusion.q)[[1]] # dmap.q(Laplacian, Markov time t, dimension q)
    Dx = as.matrix(dist((U)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc = mgc.sample(Dx, Dy, option = 'mgc')[[1]]
      mgc.stat[s] = mgc
    }
    
    if(option == 2){
      dcor = dcor.ttest(Dx, Dy, distance = TRUE)
      dcor.stat[s] = dcor$statistic
    }
      
    if(option == 3){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      hhg.stat[s] = hhg$sum.lr
    }
  }
  
  tmp = matrix(0, nrow = n.perm, ncol = length(t.range))
  
  for(r in 1:n.perm){
    
    for(s in 1:length(t.range)){
      
      diffusion.q  =  min( max(getElbows(abs(print.lambda(P, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A)-1)
      U  =  dmap.q(P, t.range[s], diffusion.q)[[1]]

      if(class(X) == "numeric"){
        per = sample(length(X));
        newX = X[per]
      }else if(class(X) == "matrix"){
        per=sample(nrow(X));
        newX = X[per,]
      }


      Dx = as.matrix(dist(U), diag = TRUE, upper = TRUE) 
      Dy = as.matrix(dist(newX), diag = TRUE, upper = TRUE)
    
      if(option == 1){
        mgc = mgc.sample(Dx, Dy, option = 'mgc')[[1]]
        tmp[r,s] = mgc
      }
      
      if(option == 2){
        dcor = dcor.ttest(Dx, Dy, distance = TRUE)
        tmp[r,s] = dcor$statistic
      }
      
      if(option == 3){
        hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
        tmp[r,s] = hhg$sum.lr
      }
    } 
  }
  
  if(option == 1) return(list(mgc.stat, tmp))
  if(option == 2) return(list(dcor.stat, tmp))
  if(option == 3) return(list(hhg.stat, tmp))
}
