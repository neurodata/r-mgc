#' Rank Data
#'
#' A function that ranks all the variables for a particular observation in a dataset.
#'
#' @param X [n, d] the data to rank, with n observations of d variables.
#' @param normalize=FALSE : a boolean indicating whether to normalize so the rank for a particular observation has values betweeen 0 and 1.
#' @return rX [n, d] the ranked data, where the ith row contains the ranks of the ith set of observations.
#' @author Eric Bridgeford
#' @export
discr.rank_data <- function(X, normalize=FALSE) {
  rX <- t(apply(X, 1, function(x) {
    # use the stats ranking function
    rx <- array(rank(x, ties.method="average"))
    if (normalize) {
      # normalize values to fall btwn 0 and 1
      rx <- ( rx - min(rx) ) / (max(rx) - min(rx))
    }
    return(rx)
  }))
  return(rX)
}
