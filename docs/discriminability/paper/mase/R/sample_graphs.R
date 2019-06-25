#' 
#' Function to sample an independent-edge graph given a probability matrix P
#' 
#' @param P matrix encoding the probabilities of an edge in the graph
#' 
#' @return The adjacency matrix of a random graph
#' 
sample_from_P <- function(P) {
  n = ncol(P)
  A = Matrix(0, n, n)
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + t(A)
  return(A)
}
