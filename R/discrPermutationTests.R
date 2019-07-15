#' Discriminability One Sample Permutation Test
#'
#' A function that permutes the labels of a distance matrix to obtain an empirical pvalue associated with whether the raw score is due to random chance.
#'
#' @importFrom parallel mclapply detectCores
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist=TRUE}.}
#' }
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples.
#' @param is.dist a boolean indicating whether your \code{X} input is a distance matrix or not. Defaults to \code{FALSE}.
#' @param dist.xfm if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{discr.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param remove.isolates remove isolated samples from the dataset. Isolated samples are samples with only
#' one instance of their class appearing in the \code{Y} vector. Defaults to \code{TRUE}.
#' @param nperm the number of permutations to perform. Defaults to \code{100}.
#' @param no_cores the number of cores to use for permutation test. Defaults to \code{1}.
#' @return A list containing the following:
#' \item{\code{srel}}{the relative, unpermuted discriminability you want to see is significant.}
#' \item{\code{null}}{the discriminability scores of the permuted data.}
#' \item{\code{pval}}{the pvalue associated with the permutation test.}
#' @author Eric Bridgeford
#' @export
discr.test.one_sample <- function(X, Y, is.dist=FALSE, dist.xfm=discr.distance, dist.params=list(method='euclidean'),
                                  dist.return=NULL, remove.isolates=TRUE, nperm=100, no_cores=1) {

  validated <- discr.validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                               remove.isolates=remove.isolates)

  D <- validated$D; Y <- validated$Y; N <- nrow(D)
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  }
  tr <- discr.stat(D, Y, is.dist=TRUE)$discr
  nr <- mclapply(1:nperm, function(i) discr.stat(D, Y[sample(N)], is.dist=TRUE)$discr,
                 no_cores=no_cores)
  result <- list()
  result$srel <- tr
  result$null <- sort(nr)
  result$pval <- (sum(nr>tr) + 1)/(nperm + 1)
  return(result)
}

#' Discriminability Two Sample Permutation Test
#'
#' A function that takes two distance matrices and produces a p-value associated with whether or not the distance matrices differ significantly.
#'
#' @importFrom parallel mclapply detectCores
#' @param X1 is interpreted as a \code{[n x d]} data matrix with \code{n} samples in \code{d} dimensions. Should NOT be a distance matrix.
#' @param X2 is interpreted as a \code{[n x d]} data matrix with \code{n} samples in \code{d} dimensions. Should NOT be a distance matrix.
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples. Should be matched such that \code{Y[i]} is the corresponding label for \code{X1[i,]} and \code{X2[i,]}.
#' @param dist.xfm if \code{is.dist == FALSE}, a distance function to transform \code{X}. If a distance function is passed,
#' it should accept an \code{[n x d]} matrix of \code{n} samples in \code{d} dimensions and return a \code{[n x n]} distance matrix
#' as the \code{$D} return argument. See \link[mgc]{discr.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param remove.isolates remove isolated samples from the dataset. Isolated samples are samples with only
#' one instance of their class appearing in the \code{Y} vector. Defaults to \code{TRUE}.
#' @param nperm the number of permutations for permutation test. Defualts to \code{100}.
#' @param no_cores the number of cores to use for the permutations. Defaults to \code{1}.
#' @param alt the alternative hypothesis. Can be that first dataset is more discriminable (\code{alt = 'greater'}), less discriminable (\code{alt = 'less'}),
#' or just non-equal (\code{alt = 'neq'}).
#' @return A list containing the following:
#' \item{\code{stat}}{the observed test statistic.}
#' \item{\code{p.value}}{The p-value associated with the test.}
#' @author Eric Bridgeford
#' @export
discr.test.two_sample <- function(X1, X2, Y, dist.xfm=discr.distance,
                                  dist.params=list(method="euclidian"), dist.return=NULL,
                                  remove.isolates=TRUE, nperm=100,
                                  no_cores=1) {

  validated1 <- discr.validator(X1, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)
  validated2 <- discr.validator(X2, Y, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                                remove.isolates=remove.isolates)
  D1 <- validated1$D; D2 <- validated2$D; Y1 <- validated1$Y; Y2 <- validated2$Y; N <- nrow(D1)
  if (!all(Y1 == Y2)) {
    stop("The ids are not equal after removal of isolates from X1 and X2.")
  }
  if (nrow(D1) != nrow(D2) || nrow(D1) != length(Y1) || nrow(D1) != ncol(D1)) {
    stop("The distance matrices do not have the same number of elements.")
  }
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  }
  # get observed D1.hat and D2.hat
  D1.hat <- discr.mnr(discr.rdf(D1, Y1)); D2.hat <- discr.mnr(discr.rdf(D2, Y1))

  null.discrs <- mclapply(1:nperm, function(i) {
    # generate null dataset for X1
    idx1 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda1 <- runif(N)
    Xn1 <- lambda1*X1[idx1[,1],] + (1 - lambda1)*X1[idx1[,2],]

    # generate null dataset for X2
    idx2 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda2 <- runif(N)
    Xn2 <- lambda2*X2[idx2[,1],] + (1 - lambda2)*X2[idx2[,2],]
    D1.null <- discr.stat(Xn1, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)
    D2.null <- discr.stat(Xn2, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)
    return(list(D1hat=D1.null, D2hat=D2.null))
  }, no_cores=no_cores)

  null.diff <- sapply(1:N, function(j) {
    sapply((j+1):N, function(j.p) {
      return(c(null.discrs[[j]]$D1hat - null.discrs[[j.p]]$D2hat, null.discrs[[j]]$D2hat - null.discrs[[j.p]]$D2hat))
    })
  })
  if (alt == 'greater') {
    stat <- D1.hat - D2.hat
    # p-value is fraction of times observed statistic is greater
    p.value <- mean(stat < null.diff)
  } else if (alt == 'less') {
    stat <- D2.hat - D1.hat
    # p-value is fraction of times observed statistic is lower
    p.value <- mean(stat < null.diff)
  } else if (alt == 'neq') {
    stat <- abs(D1.hat - D2.hat)
    # p-value is fraction of times observed statistic is more extreme
    p.value <- mean(stat < null.diff || stat > null.diff)
  } else {
    stop("You have not entered a valid alternative.")
  }
  return(list(p.value=p.value, stat=stat))
}
