#' 
#' Function to sample MASE statistic from the COSIE model. A sample of 2B graphs
#' with E[A] = VRV^T is obtained, and for every pair of them, the distance between their
#' score matrices obtained by MASE is calculated.
#' 
#' @param V common invariant subspace of COSIE
#' @param R score matrix of COSIE
#' @param P_0 optional average matrix in the model. Default is 0.
#' @param B number of test statistics samples
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NULL, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augmentation logial. When true, the diagonal augmentation method for ASE is performed.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#' 
#' @return A list containing a matrix V of size n x d, with the 
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#' 
#' @references 
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
mase.test.statistic.distr <- function(V, R, P0 = NULL, B = 1000,
                             d = NULL, d_vec = NULL, scaled.ASE = FALSE, diag.augmentation = TRUE) {
  n = nrow(V)
  d = ncol(V)
  
  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1 <- sample_from_P(P1)
    Ab2 <- sample_from_P(P1)
    mase.embedding <- mase(list(Ab1-P0, Ab2-P0), d = d, d_vec = d_vec,
                           scaled.ASE = scaled.ASE, diag.augmentation = diag.augmentation)
    return(norm(mase.embedding$R[[1]]- mase.embedding$R[[2]], "F"))
  })
  return(bootstrap)
}


#' 
#' Given a set of test statistics and a sample null distribution, calculate p-values
p.values <- function(test.statistics, null.distribution) {
  return(sapply(test.statistics, function(test) 
    (sum(test <= null.distribution) + 0.5)/ length(null.distribution)))
}


#' 
#' Function to perform a parametric bootstrap to test the hypothesis
#' that two score matrices in the COSIE model are not different.
#' 
#' @param V common invariant subspace of COSIE
#' @param R1 score matrix of COSIE for the first graph
#' @param R2 score matrix of COSIE for the first graph
#' @param B number of test statistics samples
#' @param P_0 optional average matrix in the model. Default is 0.
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NULL, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augmentation logial. When true, the diagonal augmentation method for ASE is performed.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#' 
#' @return A list containing a matrix V of size n x d, with the 
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#' 
#' @references 
#' Tang, Minh, et al. "A semiparametric two-sample hypothesis testing problem for random graphs." 
#' Journal of Computational and Graphical Statistics 26.2 (2017): 344-354.
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
test_equality_COSIE <- function(V, R1, R2, B = 1000, P0 = NULL,
                                d = NULL, d_vec = NULL, scaled.ASE = FALSE, diag.augmentation = TRUE) {
  if(is.null(P0)) P0 = Matrix(0, n, n)
  null_1 <- mase.test.statistic.distr(V, R1, P0 = P0, B = B,
                                      d = d, d_vec = d_vec, scaled.ASE = scaled.ASE, diag.augmentation = diag.augmentation)
  null_2 <- mase.test.statistic.distr(V, R2, P0 = P0, B = B,
                                      d = d, d_vec = d_vec, scaled.ASE = scaled.ASE, diag.augmentation = diag.augmentation)
  test <- norm(R1 - R2, type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2, 
              test.statistic = test))
}




