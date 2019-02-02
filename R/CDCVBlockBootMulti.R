#' Full test of independence of multivariate time series X and Y.
#' Uses cross-distance covariance and circular or stationary bootstrap.
#' For univariate data, use cdcv.univ.test.
#'
#' @param X Time series as n x d_X matrix.
#' @param Y Time series as n x d_Y matrix.
#' @param M Maximum lag to consider for cross-distance covariance.
#' @param num_boot Number of bootstrapped samples.
#' @param alpha Significance level of test.
#' @param unbiased Whether to use the biased or unbiased estimate of dCov.
#' @param type "circular" or "stationary" boostrap.
#' @param block_size Block size or mean block size.
#' @return A list containing the following:
#' \item{\code{pCDCV}}{P-value of CDCV}
#' \item{\code{statCDCV}}{test statistic value}
#' \item{\code{unbiased}}{Whether to use the unbiased or biased estimate of CDCV}
#' \item{\code{M}}{the max lead/lag considered}
#'
#'@export
cdcv.multi.test <- function(X, Y, M = NA, num_boot = 100, alpha = 0.05, unbiased = FALSE, type = "stationary", block_size = NA) {

  # Default block size is sqrt(n).
  n <- nrow(X)
  if (is.na(M)) M <- ceiling(sqrt(n))
  if (is.na(block_size)) block_size <- floor(sqrt(n))

  # Test statistic computation - a function of X holding Y fixed.
  # Must be defined internally to maintain Y in scope.
  # Uses Bartlett kernel adapted to FP statistic.
  test_stat <- function(X) {
    n <- nrow(X)
    p <- sqrt(n)
    stat_type <- if (unbiased) "U" else "V"

    result <- n*(energy::dcov(X, Y, type = stat_type))
    if (M) {
      for (j in 1:M) {
        x <- X[(1+j):n,]
        y <- Y[1:(n-j),]
        result <- result + (n-j)*((1 - j/(p*(M+1)))^2)*(energy::dcov(x, y, type = stat_type))

        x <- X[1:(n-j),]
        y <- Y[(1+j):n,]
        result <- result + (n-j)*((1 - j/(p*(M+1)))^2)*(energy::dcov(x, y, type = stat_type))
      }
    }
    result
  }

  # Block boostrap.
  if (type == "circular") sim <- "fixed" else sim <- "geom"
  boot_sample <- boot::tsboot(tseries = X,
                              statistic = test_stat,
                              R = num_boot,
                              l = block_size,
                              sim = "fixed",
                              orig.t = FALSE)

  statCDCV <- test_stat(X)
  pCDCV <- sum(as.numeric(boot_sample$t >= statCDCV))/num_boot

  result <- list(pCDCV = pCDCV,
                 statCDCV = statCDCV,
                 unbiased = unbiased,
                 max_lag = M)
  return(result)
}
