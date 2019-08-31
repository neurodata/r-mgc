#' Remove Isolates
#'
#' A function to remove isolates from a dataset, given a data matrix or a distance matrix.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x n]} distance matrix}{X is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions.}
#' }
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples.
#' @param is.dist a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @keywords internal
remove.isolates <- function(X, Y, is.dist=FALSE) {

  un.y <- as.character(unique(Y))
  id.count <- sapply(un.y, function(y) sum(Y == y))
  names(id.count) <- un.y

  for (y in un.y) {
    # if the id is an isolate, remove
    if (id.count[y] <= 1) {
      remove.rows <- which(Y == y)
      # return rows that are not the rows to be removed
      if (is.dist) {
        X <- X[-remove.rows, -remove.rows, drop=FALSE]
      } else {
        X <- X[-remove.rows,,drop=FALSE]
      }
      Y <- Y[-remove.rows]
    }
  }
  return(list(X = X, Y = Y))
}
