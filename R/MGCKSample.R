#' MGC K Sample Testing
#'
#' MGC K Sample Testing provides a wrapper for MGC Sample testing under the constraint that the Ys here are
#' categorical labels with K possible sample ids. This function uses a 0-1 loss for the Ys (one-hot-encoding)).
#'
#' @references Youjin Lee, et al. "Network Dependence Testing via Diffusion Maps and Distance-Based Correlations." ArXiv (2019).
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist.X=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist.X=TRUE}.}
#' }
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param mgc.opts Arguments to pass to MGC, as a named list. See \code{\link{mgc.test}} for details. Do not pass arguments for
#' \code{is.dist.Y}, \code{dist.xfm.Y}, \code{dist.params.Y}, nor \code{dist.return.Y}, as they will be ignored.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{p.value}}{P-value of MGC}
#' \item{\code{stat}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{pLocalCorr}}{P-value of the local correlations by double matrix index}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#' @author Eric Bridgeford
#'
#' @examples
#'
#' library(mgc)
#'
#' n = 100; d = 2
#' # simulate 200 samples which are jointly dependent in 10 dimensions
#' data <- mgc.sims.joint(n, d)
#' data_mtx <- rbind(data$X, data$Y)
#' labels <- c(replicate(n, 0), replicate(n, 1))
#' result <- mgc.ksample(data_mtx, labels, mgc.opts=list(nperm=10))
#'
#' @export
mgc.ksample <- function(X, Y, mgc.opts=list(), ...){

  DY <- mgc.distance(Y, method="ohe")

  # distance already computed, so skip...
  if (!is.null(mgc.opts$is.dist.Y)) {
    if (isFALSE(mgc.opts$is.dist.Y)) {
      warning("Warning: Distance function is automatically applied via OHE; not using user-input distance function.")
      mgc.opts$is.dist.Y <- TRUE; mgc.opts$dist.xfm.Y=NULL
      mgc.opts$dist.params.Y <- NULL; mgc.opts$dist.return.Y <- NULL
    }
  }

  return(do.call(mgc.test, c(list(X=X, Y=DY), mgc.opts)))
}
