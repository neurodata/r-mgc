#' Full test of independence of multivariate time series X and Y.
#' Uses cross-distance covariance and circular or stationary bootstrap.
#' For univariate data, use cdcv.univ.test.
#'
#' @param X Time series as n x d_X matrix.
#' @param Y Time series as n x d_Y matrix.
#' @param M Maximum lag to consider for cross-distance covariance. Defaults to \code{log(n)}.
#' @param num_boot Number of bootstrapped samples. Defaults to \code{100}.
#' @param boot_type "stationary" or "circular" bootstrap. Defaults to \code{stationary}.
#' @param block_size Block size or mean block size. Defaults to \code{sqrt(n)}.
#' @return A list containing the following:
#' \item{\code{pCDCV}}{P-value of CDCV}
#' \item{\code{statCDCV}}{test statistic value}
#' \item{\code{unbiased}}{Whether to use the unbiased or biased estimate of CDCV}
#' \item{\code{M}}{the max lead/lag considered}
#'
#'@export
cdcv.multi.test <- function(X, Y, M = NA, num_boot = 100, boot_type = "stationary", block_size = NA) {

  # Correct X and Y into matrices, if they are vectors.
  if (class(X) == "numeric") {
    X <- matrix(X, nrow = length(X))
  }
  if (class(Y) == "numeric") {
    Y <- matrix(Y, nrow = length(Y))
  }

  # Default block size is sqrt(n).
  n <- nrow(X)
  if (is.na(M)) M <- ceiling(log(n))
  if (is.na(block_size)) block_size <- floor(sqrt(n))

  # Test statistic computation - a function of X holding Y fixed.
  # Must be defined internally to maintain Y in scope.
  # Uses Bartlett kernel adapted to FP statistic.
  test_stat <- function(X) {
    p <- sqrt(n)

    result <- n*(energy::dcov(X, Y))
    if (M) {
      for (j in 1:M) {
        if (class(X) == "numeric") {
          x_lead <- matrix(X[(1+j):n], ncol = 1)
          x_lag <- matrix(X[1:(n-j)], ncol = 1)
        } else {
          x_lead <- X[(1+j):n,]
          x_lag <- X[1:(n-j),]
        }
        if (class(Y) == "numeric") {
          y_lag <- Y[1:(n-j)]
          y_lead <- Y[(1+j):n]
        } else {
          y_lag <- Y[1:(n-j),]
          y_lead <- Y[(1+j):n,]
        }

        result <- result + (n-j)*((1 - j/(p*(M+1)))^2)*(energy::dcov(x_lead, y_lag)) +
          (n-j)*((1 - j/(p*(M+1)))^2)*(energy::dcov(x_lag, y_lead))
      }
    }
    return(result)
  }

  # Block bootstrap.
  if (boot_type == "circular") sim <- "fixed" else sim <- "geom"
  boot_sample <- boot::tsboot(tseries = X,
                              statistic = test_stat,
                              R = num_boot,
                              l = block_size,
                              sim = sim,
                              orig.t = FALSE)

  statCDCV <- test_stat(X)
  pCDCV <- sum(as.numeric(boot_sample$t >= statCDCV))/num_boot

  result <- list(pCDCV = pCDCV,
                 statCDCV = statCDCV,
                 max_lag = M)
  return(result)
}
