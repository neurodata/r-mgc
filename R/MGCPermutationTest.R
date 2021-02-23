#' MGC Permutation Test
#'
#' Test of Dependence using MGC Approach.
#'
#' @references Joshua T. Vogelstein, et al. "Discovering and deciphering relationships across disparate data modalities." eLife (2019).
#' @param X is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{X is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist.X=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{X is a distance matrix. Use flag \code{is.dist.X=TRUE}.}
#' }
#' @param Y is interpreted as:
#' \describe{
#'    \item{a \code{[n x d]} data matrix}{Y is a data matrix with \code{n} samples in \code{d} dimensions, if flag \code{is.dist.Y=FALSE}.}
#'    \item{a \code{[n x n]} distance matrix}{Y is a distance matrix. Use flag \code{is.dist.Y=TRUE}.}
#' }
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
#' @param nperm specifies the number of replicates to use for the permutation test. Defaults to \code{1000}.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{use the MGC global correlation.}
#'    \item{\code{'dcor'}}{use the dcor global correlation.}
#'    \item{\code{'mantel'}}{use the mantel global correlation.}
#'    \item{\code{'rank'}}{use the rank global correlation.}
#' }
#' @param no_cores the number of cores to use for the permutations. Defaults to \code{1}.
#'
#' @return A list containing the following:
#' \item{\code{p.value}}{P-value of MGC}
#' \item{\code{stat}}{is the sample MGC statistic within \code{[-1,1]}}
#' \item{\code{p.localCorr}}{P-value of the local correlations by double matrix index.}
#' \item{\code{localCorr}}{the local correlations}
#' \item{\code{optimalScale}}{the optimal scale identified by MGC}
#' \item{\code{option}}{specifies which global correlation was used}
#'
#' @author Eric Bridgeford and C. Shen
#'
#' @section Details:
#'
#' A test of independence using the MGC approach, described in Vogelstein et al. (2019). For \eqn{X \sim F_X}{X ~ Fx}, \eqn{Y \sim F_Y}{Y ~ Fy}:
#'
#' \deqn{H_0: F_X \neq F_Y}{H0: Fx != Fy} and: \deqn{H_A: F_X = F_Y}{Ha: Fx = Fy}
#'
#' Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
#' unless corrected for multiple hypotheses.
#'
#' For details on usage see the help vignette:
#' \code{vignette("mgc", package = "mgc")}
#'
#'
#' @examples
#' \dontrun{
#' library(mgc)
#'
#' n = 100; d = 2
#' data <- mgc.sims.linear(n, d)
#' # note: on real data, one would put nperm much higher (at least 100)
#' # nperm is set to 10 merely for demonstration purposes
#' result <- mgc.test(data$X, data$Y, nperm=10)
#'}
#' @export
mgc.test <-function(X, Y, is.dist.X=FALSE, dist.xfm.X=mgc.distance, dist.params.X=list(method='euclidean'),
                    dist.return.X=NULL, is.dist.Y=FALSE, dist.xfm.Y=mgc.distance, dist.params.Y=list(method='euclidean'),
                    dist.return.Y=NULL, nperm=1000, option='mgc', no_cores=1) {
  # validate input is valid and convert to distance matrices, if necessary
  validated <- mgc.validator(X, Y, is.dist.X=is.dist.X, dist.xfm.X=dist.xfm.X, dist.params.X=dist.params.X,
                             dist.return.X=dist.return.X, is.dist.Y=is.dist.Y, dist.xfm.Y=dist.xfm.Y, dist.params.Y=dist.params.Y,
                             dist.return.Y=dist.return.Y)

  DX <- validated$DX; DY <- validated$DY

  if (no_cores > detectCores()) {
    stop(sprintf("Requested more cores than available. Requested %d cores; CPU has %d.", no_cores, detectCores()))
  } else if (no_cores >= 0.8*detectCores()) {
    cat("You have requested a number of cores near your machine's core count. Expected decreased performance.\n")
  }
  N = nrow(DY)

  if (nperm < 100) {
    warning("nperm is < 100. nperm should typically be set > 100.")
  }

  if (!min(abs(c(nperm%%1, nperm%%1-1))) < 1e-5) {
    stop("You have not specified an integer for `nperm`.")
  }

  # Compute sample MGC and all local correlations
  result <- mgc.stat.driver(DX, DY, option)

  # compute the null using permutation test
  mgc.nulls <- mclapply(1:nperm, function(i) {
    per <- sample(N)  # resample the IDs for Y
    return(mgc.stat.driver(DX, DY[per, per], option=option))
  }, mc.cores=no_cores)

  # get the boolean mgc test results for each permutation
  result <- list(stat=result$stat, localCorr=result$localCorr, optimalScale=result$optimalScale,
                 option=option)

  #  if (isFALSE(fast)) {
  pLocalCorrs <- lapply(mgc.nulls, function(mgc.null) mgc.null$localCorr >= result$localCorr)
  result$p.localCorr <- (Reduce("+", pLocalCorrs) + 1)/(nperm + 1)
  #
  pMGCs <- lapply(mgc.nulls, function(mgc.null) mgc.null$stat >= result$stat)

  # element-wise average
  result$p.value <- (Reduce("+", pMGCs) + 1)/(nperm + 1)

  return(result)
}

