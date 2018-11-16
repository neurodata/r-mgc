#' Reliability Density Function
#'
#' A function for computing the reliability density function of a dataset.
#'
#' @param X \code{[n, n]} a distance matrix for n samples.
#' @param ids \code{[n]} a vector containing the label ids for each sample.
#' @return \code{[n]} vector of the reliability per sample.
#' @author Eric Bridgeford
discr.rdf <- function(X, ids) {
  N <- length(ids)  # number of scans total
  K <- unique(ids)  # unique subjects
  # some basic checks
  # valid valid N
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  # validate X is a distance matrix
  if (!all.equal(dim(X), c(N, N))) {
    stop('You have not passed a valid datatype for X.')
  }

  # loop over scans
  rdf <- sapply(1:N, function(i) {
    ind <- which(ids[i] == ids) # all the indices that are the same subject, but different scan
    di <- X[i,] # get the entire ra for the particular scan
    # loop over the unique subjects
    inn <- sapply(K, function(j) {
      indj <- which(j == ids)  # get index of subject j's scans
      Dij <- di[indj[indj != i]]  # indices of D corresponding to between scan i and subject j, excluding scan i
      return(min(Dij, na.rm=TRUE))  # nearest neighbor from scan i to collection of scans for subject j
    })
    # d between current scan and nearest neighbord of current subject
    d <- inn[K == ids[i]]
    # ds between current scan and nearest neighbor of other subject
    inn <- inn[K != ids[i]]
    return(1 - (sum(inn[!is.nan(inn)] < d) + 0.5*sum(inn[!is.nan(inn)] == d)) / (length(K) - 1))
  })

  return(rdf)
}

#' Discriminability Statistic
#'
#' A function for computing the discriminability from a distance matrix and a set of associated labels.
#'
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x n]} distance matrix}{X is a square matrix with zeros on diagonal for \code{n} samples.}
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions.}
#' }
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples.
#' @param remove_outliers boolean indicating whether to ignore observations with rdf below a certain cutoff. Defaults to \code{FALSE}.
#' @param thresh the threshold below which to ignore observations. If thresh > 0, ignores observations where the rdf is < thresh in the discriminability computation. Defaults to \code{0}.
#' @param verbose a boolean indicating whether to print output. Defaults to \code{FALSE}.
#' @return A list containing the following:
#' \item{`discr` the discriminability statistic.}
#' \item{`rdf` the rdfs for each sample.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("discriminability", package = "mgc")}
#'
#' @examples
#'
#' nsrc <- 5
#' nobs <- 10
#' d <- 20
#' set.seed(12345)
#' src_id <- array(1:nsrc)
#' labs <- sample(rep(src_id, nobs))
#' dat <- t(sapply(labs, function(lab) rnorm(d, mean=lab, sd=1)))
#' discr.stat(dat, labs)
#'
#' @author Eric Bridgeford
#' @export
discr.stat <- function(X, Y) {
  X <- as.matrix(X)
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(X) != ncol(X) | sum(diag(X)^2) > 0){
    X <- discr.distance(X)
  }
  rdf <- discr.rdf(X, Y)
  return(list(discr=mean(rdf, is.nan=FALSE), rdf=rdf))
}
