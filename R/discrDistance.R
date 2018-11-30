#' Distance
#'
#' A function that returns a distance matrix given a collection of observations.
#'
#' @importFrom stats dist
#' @param X \code{[n x d]} a data matrix for \code{d} samples of \code{d} variables.
#' @param method the method for computing distances. Defaults to \code{'euclidean'}. See \link[stats]{dist} for details.
#' @return a \code{[n x n]} distance matrix.
#' \code{[n, n]} distance matrix indicating the pairwise distances between all samples passed in.
#' @author Eric Bridgeford
#' @export
discr.distance <- function(X, method='euclidean') {
  D <- as.matrix(dist(X, diag=TRUE, method=method))
  return(D)
}
