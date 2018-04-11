#' Distance
#'
#' A function that returns a distance matrix given a collection of observations.
#'
#' @importFrom stats dist
#' @param X \code{[n, d]} a data matrix for n samples of d variables.
#' @param method the method for computing distances. Defaults to \code{'2'}.
#' \describe{
#'    \item{'2'}{use the 2 norm between pairs of observations}
#' }
#' @return a \code{[n, n]} distance matrix indicating the pairwise distances between all samples passed in.
#' @author Eric Bridgeford
discr.distance <- function(X, method='2') {
  if (method == '2') {
    D <- as.matrix(dist(X, diag=TRUE))
  } else {
    stop('You have not passed a valid method. Please select one of the following: [\'2\'].')
  }
  return(D)
}
