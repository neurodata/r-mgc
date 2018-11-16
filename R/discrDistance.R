#' Distance
#'
#' A function that returns a distance matrix given a collection of observations.
#'
#' @importFrom stats dist
#' @param X \code{[n, d]} a data matrix for `n` samples of `d` variables.
#' @param method the method for computing distances. Defaults to \code{'euclidian'}. See `\link[stats]{dist}` for details.
#' @return a \code{[n, n]} distance matrix indicating the pairwise distances between all samples passed in.
#' @author Eric Bridgeford
discr.distance <- function(X, method='euclidian') {
  if (!is.matrix(X)) {
    tryCatch({
      X <- as.matrix(X)
    }, error=function(e) stop("You have not passed an object that can be coerced to a '[n, d] matrix'."))
  }
  D <- as.matrix(dist(X, diag=TRUE, method=method))

  return(D)
}
