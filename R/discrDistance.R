#' Distance
#'
#' A function that returns a distance matrix given a collection of observations.
#'
#' @import stats
#' @param X [n, d] a data matrix for n observations of d variables.
#' @param method='2' the method for computing distances.
#' \describe{
#'    \item{'2'}{use the 2 norm between pairs of observations}
#' }
#' @return dist [n, n]: a matrix indicating the pairwise distances between all observations passed in.
#' @author Eric Bridgeford
#' @export
discr.distance <- function(X, method='2') {
  if (method == '2') {
    D <- as.matrix(dist(X, diag=TRUE))
  } else {
    stop('You have not passed a valid method. Please select one of the following: [\'2\'].')
  }
  return(D)
}
