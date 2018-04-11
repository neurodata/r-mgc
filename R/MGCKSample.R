#' MGC K Sample Testing
#'
#' MGC K Sample Testing provides a wrapper for MGC Sample testing under the constraint that the Ys here are
#' categorical labels with K possible sample ids. This function uses a 0-1 loss for the Ys. To use a custom distance
#' function, use \code{mgc.test} with your custom distance function as \code{Y}.
#'
#' @param X is interpreted as:
#' \itemize{
#'    \item{\code{[n x n]} distance matrix}{X is a square matrix with zeros on diagonal}
#'    \item{\code{[n x d]} data matrix}{Otherwise}
#'  }
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param mgc.opts Arguments to pass to MGC. See \code{\link{mgc.test}} for details.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{pMGC}}{P-value of MGC}
#' \item{\code{statMGC}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{pLocalCorr}}{P-value of the local correlations by double matrix index}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#' @author Eric Bridgeford
#'
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
#' result <- mgc.ksample(data_mtx, labels, mgc.opts=list(rep=10))
#'
#' @export
mgc.ksample = function(X, Y, mgc.opts=list(), ...){
  ylabs <- unique(Y); K <- length(ylabs); n <- length(Y)
  # one-hot-encode the y-labels for 0-1 loss under euclidian distance
  Yh <- array(0, dim=c(n, K))
  for (i in 1:K) {
    Yh[Y == ylabs[i],i] <- 1
  }

  # compute distance...
  Dy <- as.matrix(dist(Yh, method='binary'), nrow=n)

  result =  do.call(mgc.test, c(list(X=X, Y=Dy), mgc.opts))

  return(result)
}
