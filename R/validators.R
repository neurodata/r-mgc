#' Discriminability Utility Validator
#'
#' A script that validates that data inputs are correct, and returns a distance matrix and a ids vector.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param Y is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{Y is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{Y is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param is.dist a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param remove.isolates whether to remove isolated samples, or samples with only a single instance in the \code{Y} vector.
#' @return A list containing the following:
#' \item{\code{DX}}{The X distance matrix, as a \code{[n x n]} matrix.}
#' \item{\code{Y}}{The sample ids, with isolates removed.}
#' @keywords internal
discr.validator <- function(X, Y, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                            dist.return=NULL, remove.isolates=TRUE) {

  DX <- mgc.dist.validator(X, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return)

  # validate Y
  if (!is.vector(Y)) {
    tryCatch({
      Y <- as.vector(Y)
    }, error=function(e) stop("You have not passed an object Y that can be coerced to a [n] vector."))
  }

  if (nrow(DX) != length(Y)) {
    stop("Your distance matrix and your ids vector do not have the same number of samples, after applying distance function.")
  }

  # remove isolated subjects if requested.
  if (remove.isolates) {
    purged <- remove.isolates(DX, Y, is.dist=TRUE)
    DX <- purged$X; Y <- purged$Y
  }


  if (length(unique(Y)) <= 1) {
    stop("You have passed a vector containing only a single unique sample id.")
  }

  return(list(D=DX, Y=Y))
}

#' MGC Utility Validator
#'
#' A script that validates that data inputs are correct, and returns a X distance and Y distance matrix for MGC.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples.
#' @param is.dist.X a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm.X if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params.X a list of trailing arguments to pass to the distance function specified in \code{dist.xfm.X}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return.X the return argument for the specified \code{dist.xfm.X} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm.X[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param is.dist.Y a boolean indicating whether your \code{Y} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm.Y if \code{is.dist == FALSE}, a distance function to transform \code{Y}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{dist.return.Y} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params.Y a list of trailing arguments to pass to the distance function specified in \code{dist.xfm.Y}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return.Y the return argument for the specified \code{dist.xfm.Y} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm.Y(Y)} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm.Y(Y)[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @return A list containing the following:
#' \item{\code{D}}{The distance matrix, as a \code{[n x n]} matrix.}
#' \item{\code{Y}}{the sample ids, as a \code{[n]} vector.}
#' @keywords internal
mgc.validator <- function(X, Y, is.dist.X=FALSE, dist.xfm.X=mgc.distance, dist.params.X=list(method='euclidean'),
                          dist.return.X=NULL, is.dist.Y=FALSE, dist.xfm.Y=mgc.distance, dist.params.Y=list(method='euclidean'),
                          dist.return.Y=NULL) {

  DX <- mgc.dist.validator(X, is.dist=is.dist.X, dist.xfm=dist.xfm.X, dist.params=dist.params.X, dist.return=dist.return.X)
  DY <- mgc.dist.validator(Y, is.dist=is.dist.Y, dist.xfm=dist.xfm.Y, dist.params=dist.params.Y, dist.return=dist.return.Y)

  if (nrow(X) != nrow(Y)) {
    stop(sprintf("Your X and Y distance functions produced distance matrices with different #samples. X has %d samples; Y has %d samples.",
                 dim(X)[1], dim(Y)[1]))
  }

  return(list(DX=DX, DY=DY))
}

#' Distance Matrix Validator
#'
#' A utility to validate a distance matrix.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param is.dist a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @return A distance matrix.
#' @author Eric Bridgeford
#' @keywords internal
mgc.dist.validator <- function(X, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                               dist.return=NULL) {
  if (is.dist) {
    if (nrow(X) != ncol(X)) {
      stop("Your data X is not a [n x n] distance matrix.")
    }
  }

  if (!is.matrix(X)) {
    tryCatch({
      X <- as.matrix(X)
    }, error=function(e) stop("You have not passed an object that can be coerced to a [n x d] matrix."))
  }

  # Distance transform if requested by the user.
  if (!is.dist) {
    tryCatch({
      X <- do.call(dist.xfm, c(list(X), dist.params))
      if (!is.null(dist.return)) {
        X <- X[[dist.return]]
      }
    }, error=function(e) {
      print("Your distance function requested experienced an error.")
      stop(e)
    })
    if (nrow(X) != ncol(X)) {
      stop(sprintf("Your distance function returned an invalid distance matrix. The return should be [n x n]; yours is [%d x %d].",
                   dim(X)[1], dim(X)[2]))
    }
  }

  return(X)
}
