#' Discriminability One Sample Permutation Test
#'
#' A function that performs a one-sample test for whether the discriminability differs from random chance.
#'
#' @references Eric W. Bridgeford, et al. "Optimal Decisions for Reference Pipelines and Datasets: Applications in Connectomics." Bioarxiv (2019).
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
#' as the \code{$D} return argument. See \link[mgc]{mgc.distance} for details.
#' @param dist.params a list of trailing arguments to pass to the distance function specified in \code{dist.xfm}.
#' Defaults to \code{list(method='euclidean')}.
#' @param dist.return the return argument for the specified \code{dist.xfm} containing the distance matrix. Defaults to \code{FALSE}.
#' \describe{
#'     \item{\code{is.null(dist.return)}}{use the return argument directly from \code{dist.xfm} as the distance matrix. Should be a \code{[n x n]} matrix.}
#'     \item{\code{is.character(dist.return) | is.integer(dist.return)}}{use \code{dist.xfm[[dist.return]]} as the distance matrix. Should be a \code{[n x n]} matrix.}
#' }
#' @param remove.isolates remove isolated samples from the dataset. Isolated samples are samples with only
#' one instance of their class appearing in the \code{Y} vector. Defaults to \code{TRUE}.
#' @param nperm the number of permutations to perform. Defaults to \code{500}.
#' @param no_cores the number of cores to use for permutation test. Defaults to \code{1}.
#' @return A list containing the following:
#' \item{\code{stat}}{the discriminability of the data.}
#' \item{\code{null}}{the discriminability scores under the null, computed via permutation.}
#' \item{\code{p.value}}{the pvalue associated with the permutation test.}
#' @author Eric Bridgeford
#'
#' @section Details:
#'
#' Performs a test of whether an observed discriminability is significantly different from chance, as described in Bridgeford et al. (2019).
#' With \eqn{\hat D_X}{Dhatx} the sample discriminability of \eqn{X}{X}:
#' \deqn{H_0: D_X = D_0}{H0: Dx = D0} and:\deqn{H_A: D_X > D_0}{Ha: Dx > D0} where \eqn{D_0}{D0}
#' is the discriminability that would be observed by random chance.
#'
#'
#' @examples
#' \dontrun{
#' require(mgc)
#' n = 100; d=5
#'
#' # simulation with a large difference between the classes
#' # meaning they are more discriminable
#' sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=10)
#' X <- sim$X; Y <- sim$Y
#'
#' # p-value is small
#' discr.test.one_sample(X, Y)$p.value
#'}
#' @export
discr.test.one_sample <- function(X, Y, is.dist=FALSE, dist.xfm=mgc.distance, dist.params=list(method='euclidean'),
                                  dist.return=NULL, remove.isolates=TRUE, nperm=500, no_cores=1) {

  validated <- discr.validator(X, Y, is.dist=is.dist, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                               remove.isolates=remove.isolates)

  D <- validated$D; Y <- validated$Y; N <- nrow(D)
  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  } else if (no_cores >= 0.8*detectCores()) {
    warning("You have requested a number of cores near your machine's core count. Expected decreased performance.")
  }
  tr <- discr.stat(D, Y, is.dist=TRUE)$discr
  nr <- unlist(mclapply(1:nperm, function(i) {
    return(discr.stat(D, Y[sample(N)], is.dist=TRUE)$discr)
    }, mc.cores=no_cores), use.names=FALSE)
  result <- list()
  result$stat <- tr
  result$null <- sort(nr)
  result$p.value <- (sum(nr>tr) + 1)/(nperm + 1)
  return(result)
}

#' Discriminability Two Sample Permutation Test
#'
#' A function that takes two sets of paired data and tests of whether or not the data is more, less, or non-equally discriminable between the set of paired data.
#'
#' @references Eric W. Bridgeford, et al. "Optimal Decisions for Reference Pipelines and Datasets: Applications in Connectomics." Bioarxiv (2019).
#' @importFrom parallel mclapply detectCores
#' @param X1 is interpreted as a \code{[n x d]} data matrix with \code{n} samples in \code{d} dimensions. Should NOT be a distance matrix.
#' @param X2 is interpreted as a \code{[n x d]} data matrix with \code{n} samples in \code{d} dimensions. Should NOT be a distance matrix.
#' @param Y \code{[n]} a vector containing the sample ids for our \code{n} samples. Should be matched such that \code{Y[i]} is the corresponding label for \code{X1[i,]} and \code{X2[i,]}.
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
#' @param remove.isolates remove isolated samples from the dataset. Isolated samples are samples with only
#' one instance of their class appearing in the \code{Y} vector. Defaults to \code{TRUE}.
#' @param nperm the number of permutations for permutation test. Defualts to \code{500}.
#' @param no_cores the number of cores to use for the permutations. Defaults to \code{1}.
#' @param alt the alternative hypothesis. Can be that first dataset is more discriminable (\code{alt = 'greater'}), less discriminable (\code{alt = 'less'}),
#' or just non-equal (\code{alt = 'neq'}). Defaults to \code{"greater"}.
#' @return A list containing the following:
#' \item{\code{stat}}{the observed test statistic. the test statistic is the difference in discriminability of X1 vs X2.}
#' \item{\code{discr}}{the discriminabilities for each of the two data sets, as a list.}
#' \item{\code{null}}{the null distribution of the test statistic, computed via permutation.}
#' \item{\code{p.value}}{The p-value associated with the test.}
#' \item{\code{alt}}{The alternative hypothesis for the test.}
#' @author Eric Bridgeford
#'
#' @section Details:
#'
#' A function that performs a two-sample test for whether the discriminability is different for that of
#' one dataset vs another, as described in Bridgeford et al. (2019). With \eqn{\hat D_{X_1}}{Dhatx1} the sample discriminability of one approach, and \eqn{\hat D_{X_2}}{Dhatx2} the sample discriminability of another approach:
#'
#' \deqn{H_0: D_{X_1} = D_{X_2}}{H0: Dx1 = Dx2} and:\deqn{H_A: D_{X_1} > D_{X_2}}{Ha: Dx1 > Dx2}.
#' Also implemented are tests of \eqn{<}{<} and \eqn{\neq}{!=}.
#'
#'
#' @examples
#' \dontrun{
#' require(mgc)
#' require(MASS)
#'
#' n = 100; d=5
#'
#' # generate two subjects truths; true difference btwn
#' # subject 1 (column 1) and subject 2 (column 2)
#' mus <- cbind(c(0, 0), c(1, 1))
#' Sigma <- diag(2)  # dimensions are independent
#'
#' # first dataset X1 contains less noise than X2
#' X1 <- do.call(rbind, lapply(1:dim(mus)[2],
#'   function(k) {mvrnorm(n=50, mus[,k], 0.5*Sigma)}))
#' X2 <- do.call(rbind, lapply(1:dim(mus)[2],
#'   function(k) {mvrnorm(n=50, mus[,k], 2*Sigma)}))
#' Y <- do.call(c, lapply(1:2, function(i) rep(i, 50)))
#'
#' # X1 should be more discriminable, as less noise
#' discr.test.two_sample(X1, X2, Y, alt="greater")$p.value  # p-value is small
#' }
#' @export
discr.test.two_sample <- function(X1, X2, Y, dist.xfm=mgc.distance,
                                  dist.params=list(method="euclidian"), dist.return=NULL,
                                  remove.isolates=TRUE, nperm=500,
                                  no_cores=1, alt="greater") {

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
  } else if (no_cores >= 0.8*detectCores()) {
    warning("You have requested a number of cores near your machine's core count. Expected decreased performance.")
  }
  # get observed D1.hat and D2.hat
  D1.hat <- discr.mnr(discr.rdf(D1, Y1)); D2.hat <- discr.mnr(discr.rdf(D2, Y1))

  null.discrs <- mclapply(1:nperm, function(i) {
    # generate null dataset for X1
    idx1 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda1 <- runif(N)
    Xn1 <- lambda1*X1[idx1[,1],] + (1 - lambda1)*X1[idx1[,2],]  # convex combination of elements of X1

    # generate null dataset for X2
    idx2 <- t(sapply(1:N, function(j) sample(N, size=2)))
    lambda2 <- runif(N)
    Xn2 <- lambda2*X2[idx2[,1],] + (1 - lambda2)*X2[idx2[,2],]  # convex combination of elements of X2

    # compute discriminability under the null
    D1.null <- discr.stat(Xn1, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr
    D2.null <- discr.stat(Xn2, Y1, is.dist=FALSE, dist.xfm=dist.xfm, dist.params=dist.params, dist.return=dist.return,
                          remove.isolates=remove.isolates)$discr
    return(list(D1hat=D1.null, D2hat=D2.null))
  }, mc.cores=no_cores)

  # compute null distribution of difference between discriminabilities
  null.diff <- do.call(c, sapply(1:(nperm-1), function(j) {
    as.vector(sapply((j+1):nperm, function(j.p) {
      return(c(null.discrs[[j]]$D1hat - null.discrs[[j.p]]$D2hat, null.discrs[[j]]$D2hat - null.discrs[[j.p]]$D1hat))
    }))
  }))
  stat <- D1.hat - D2.hat
  if (alt == 'greater') {
    # p-value is fraction of times observed statistic is less
    # than under null
    p.value <- mean(stat < null.diff)
  } else if (alt == 'less') {
    # p-value is fraction of times observed statistic is greater
    # than under null
    p.value <- mean(stat > null.diff)
  } else if (alt == 'neq') {
    # p-value is fraction of times observed statistic is less extreme
    # than under null
    p.value <- mean(abs(stat) < abs(null.diff))
  } else {
    stop("You have not entered a valid alternative.")
  }
  return(list(p.value=p.value*(nperm)/(nperm + 1) + 1/(nperm + 1), stat=stat,
              discr=list(X1=D1.hat, X2=D2.hat), null=null.diff, alt=alt))
}
