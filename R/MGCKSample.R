#' MGC K Sample Testing
#'
#' MGC K Sample Testing provides a wrapper for MGC Sample testing under the constraint that the Ys here are
#' categorical labels with K possible sample ids. This function uses a 0-1 loss for the Ys. To use a custom distance
#' function, use \code{mgc.test} with your custom distance function as \code{Y}.
#'
#' @param X \code{[n, d]} the data with \code{n} samples in \code{d} dimensions.
#' @param Y \code{[n]} the labels of the samples with \code{K} unique labels.
#' @param mgc.opts Arguments to pass to MGC. See \code{\link{MGC}} for details.
#' @param ... trailing args.
#' @return A list containing the following:
#' \item{\code{pMGC}}{P-value of MGC}
#' \item{\code{statMGC}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{pLocalCorr}}{P-value of the local correlations by double matrix index}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#' @author Eric Bridgeford
#' @export
mgc.ksample = function(X, Y, mgc.opts=list(), ...){
  ylabs <- unique(Y); K <- length(ylabs)
  # one-hot-encode the y-labels for 0-1 loss under euclidian distance
  Yh <- array(0, dim=c(n, K))
  for (i in 1:K) {
    Yh[Y == ylabs[i],i] <- 1
  }

  # compute distance...
  Dy <- dist(Yh)

  result =  do.call(mgc.test, c(list(X=Dx, Y=Dy), mgc.opts))

  return(result)
}
