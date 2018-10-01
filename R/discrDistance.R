#' Distance
#'
#' A function that returns a distance matrix given a collection of observations.
#'
#' @importFrom stats dist
#' @param X \code{[n, d]} a data matrix for n samples of d variables.
#' @param method the method for computing distances. Defaults to \code{'euclidian'}.
#' \describe{
#'    \item{'2'}{use the 2 norm between pairs of observations}
#' }
#' @return a \code{[n, n]} distance matrix indicating the pairwise distances between all samples passed in.
#' @author Eric Bridgeford
discr.distance <- function(X, method='euclidian') {
  D <- as.matrix(dist(X, diag=TRUE, method=method))

  return(D)
}
