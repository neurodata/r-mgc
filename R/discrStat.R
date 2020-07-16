#' Discriminability Statistic
#'
#' A function for computing the discriminability from a distance matrix and a set of associated labels.
#'
#' @references Eric W. Bridgeford, et al. "Optimal Decisions for Reference Pipelines and Datasets: Applications in Connectomics." Bioarxiv (2019).
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples.
#' @param is.dist a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix,
#' which can be either the default output, an item castable to a distance matrix, or . See \link[mgc]{mgc.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be an object castable to a \code{[n x n]} matrix.
#'     You can verify whether this is the case by looking at \code{as.matrix(do.call(dist.xfm, list(X, <trailing_args>))}}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param remove.isolates remove isolated samples from the dataset. Isolated samples are samples with only
#' one instance of their class appearing in the \code{Y} vector. Defaults to \code{TRUE}.
#' @return A list containing the following:
#' \item{\code{discr}}{the discriminability statistic.}
#' \item{\code{rdf}}{the rdfs for each sample.}
#'
#' @section Details:
#' For more details see the help vignette:
#' \code{vignette("discriminability", package = "mgc")}
#'
#' @examples
#' sim <- discr.sims.linear(100, 10, K=2)
#' X <- sim$X; Y <- sim$Y
#' discr.stat(X, Y)$discr
#'
#' @author Eric Bridgeford
#' @export
discr.stat <- function(X, Y, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                       dist.return=NULL, remove.isolates=TRUE) {
  validated <- discr.validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                               remove.isolates=remove.isolates)
  D <- validated$D; Y <- validated$Y
  rdf <- discr.rdf(D, Y)
  return(list(discr=discr.mnr(rdf), rdf=rdf))
}


#' Reliability Density Function
#'
#' A function for computing the reliability density function of a dataset.
#'
#' @param X \code{[n, n]} a distance matrix for n samples.
#' @param ids \code{[n]} a vector containing the label ids for each sample.
#' @return \code{[n]} vector of the reliability per sample.
#' @keywords internal
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
  rdf <- sapply(1:length(ids), function(i) {
    ind <- which(ids[i] == ids) # all the indices that are the same subject
    di <- X[i,]
    Dii <- di[ind[ind != i]]  # indices of D corresponding to between scan i subject of scan i, excluding scan i
    Dij <- di[-c(ind)]  # indices of D corresponding to between i and those of other subjects, excluding scans of this subject
    # if there exist scans associated with this subject othe rthan current scan and
    # other subjects in dataset
    if (length(Dii) > 0 & length(Dij) > 0) {
      # return discriminability local rdfs as an array
      return(unlist(sapply(Dii, function(d) {
        # discriminability is (1 - count(< d) + 0.5*(count == d))/(N) where N is number of samples with different sample ids
        1 - (sum(as.numeric(Dij < d)) + 0.5*sum(Dij == d))/length(Dij)
      })))
    } else {
      warning(sprintf("Id %s is isolated in your dataset. Consider setting 'remove.isolates = TRUE'.", ids[i]))
      return(NaN)
    }
  })

  return(unlist(rdf))
}

#' Discriminability Mean Normalized Rank
#' @param rdf the reliability densities.
#' @keywords internal
#' @return the mnr.
discr.mnr <- function(rdf) {
  mean(rdf, is.nan=FALSE)
}
