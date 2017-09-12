#' Hellinger Distances
#'
#' function that computes the hellinger distance between two matrices
#'
#' @author Eric Bridgeford
#' @keywords hellinger distance
#' @param a the first matrix
#' @param b the second matrix
#' @return h the hellinger distance
#' @export
discr.hell_dist <- function(a, b) {
  return(1/sqrt(2)*norm(sqrt(a) - sqrt(b), "f"))
}
