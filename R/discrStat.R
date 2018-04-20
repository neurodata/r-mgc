#' Reliability Density Function
#'
#' A function for computing the reliability density function of a dataset.
#'
#' @param D \code{[n, n]} a distance matrix for n samples.
#' @param ids \code{[n]} a vector containing the label ids for each sample.
#' @return \code{[n]} vector of the reliability per sample.
#' @author Eric Bridgeford
discr.rdf <- function(D, ids) {
  N <- dim(D)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }

  uniqids <- unique(as.character(ids))
  countvec <- vector(mode="numeric",length=length(uniqids))

  for (i in 1:length(uniqids)) {
    countvec[i] <- sum(uniqids[i] == ids) # total number of scans for the particular id
  }

  scans <- max(countvec) # assume that we will worst case have the most
  rdf <- array(NaN, N*(scans-1)) # initialize empty ra

  count <- 1
  for (i in 1:N) {
    ind <- which(ids[i] == ids) # all the indices that are the same subject, but different scan
    for (j in ind) {
      if (!isTRUE(all.equal(j, i))) { # if j != i, then we want j to have a close distance to i, and estimate where it ranks
        di <- D[i,] # get the entire ra for the particular scan
        di[ind] <- Inf # don't want to consider the particular scan itself
        d <- D[i,j] # the distance between the particular scan of a subject and another scan of the subject
        rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind)) # 1 for less than, .5 if equal, then average
        count <-  count + 1
      }
    }
  }
  return(rdf[1:count-1]) # return only the occupied portion
}

#' Mean Normalized Rank
#'
#' A function for computing the mnr from an rdf.
#'
#' @param rdf \code{[n]} the reliability density function for \code{n} samples.
#' @param remove_outliers boolean indicating whether to ignore samples with rdf below a certain cutoff. Defaults to \code{FALSE}.
#' @param thresh the threshold below for \code{rdf} which to ignore samples. Defaults to \code{0}.
#' @param output a boolean indicating whether to ignore output. Defaults to \code{False}.
#' @return \code{discr} the discriminability statistic.
#' @author Eric Bridgeford
discr.mnr <- function(rdf, remove_outliers=FALSE, thresh=0, output=FALSE) {
  if (remove_outliers) {
    discr <- mean(rdf[which(rdf[!is.nan(rdf)] > thresh)]) # mean of the rdf
    ol <- length(which(rdf<thresh))
    if (output) {
      print(paste('Graphs with reliability <',thresh,'(outliers):', ol))
    }
  } else {
    ol <- 0
    discr <- mean(rdf[!is.nan(rdf)])
  }
  nopair <- length(rdf[is.nan(rdf)])
  if (output) {
    print(paste('Graphs with unique ids:',nopair))
    print(paste('Graphs available for reliability analysis:', length(rdf)-ol-nopair))
    print(paste('discr:', discr))
  }
  return(discr)
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
#' @param ids \code{[n]} a vector containing the labels for our \code{n} samples.
#' @param remove_outliers boolean indicating whether to ignore observations with rdf below a certain cutoff. Defaults to \code{FALSE}.
#' @param thresh the threshold below which to ignore observations. If thresh > 0, ignores observations where the rdf is < thresh in the discriminability computation. Defaults to \code{0}.
#' @param verbose a boolean indicating whether to print output. Defaults to \code{FALSE}.
#' @return discr the discriminability statistic.
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
discr.stat <- function(X, ids, remove_outliers=FALSE, thresh=0, verbose=FALSE) {
  X <- as.matrix(X)
  # Use the data size and diagonal element to determine if the given data is a distance matrix or not
  if (nrow(X) != ncol(X) | sum(diag(X)^2) > 0){
    X <- discr.distance(X)
  }
  return(discr.mnr(discr.rdf(X, ids), remove_outliers=remove_outliers, thresh=thresh, output=(verbose)))
}
