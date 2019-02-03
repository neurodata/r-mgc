context("cdcv-univ")

#' Bivariate VARMA process simulation.
#'
#' @param n Sample size.
#' @param Phis A list of 2 by 2 matrices, representing Phi coeffs.
#' @param Thetas A list of 2 by 2 matrices, representing Theta coeffs.
#' @param Sigma 2 x 2 covariance matrix for noise.
#' @return n x 2 matrix of observations, x dimensions then y dimensions.
varma <- function(n, Phis, Thetas, Sigma) {
  d <- nrow(Sigma)
  Z <- matrix(rep(0, d*n), n, d)
  innov <- MASS::mvrnorm(n, rep(0, d), Sigma)
  Z[1,] <- innov[1,]
  for (t in 2:n) {
    # AR.
    if (length(Phis)) {
      for (i in 1:length(Phis)) {
        if (t-i > 0) {
          Z[t,] <- Z[t,] + Z[t-i,] %*% Phis[[i]]
        }
      }
    }
    # MA.
    if (length(Thetas)) {
      for (j in 1:length(Thetas)) {
        if (t-j > 0) {
          Z[t,] <- Z[t,] + innov[t-j,] %*% Thetas[[j]]
        }
      }
    }
    # Innovation.
    Z[t,] <- Z[t,] + innov[t,]
  }
  colnames(Z) <- c("X", "Y")
  Z
}

test_that("All 0 inputs.", {
  n <- 20
  X <- rep(0, n)
  Y <- rep(0, n)
  result <- cdcv.univ.test(X,Y)
  expect_equal(result$pCDCV, 1)
  expect_equal(result$statCDCV, 0)
  expect_false(result$unbiased)
  expect_equal(result$max_lag, 5)
})

test_that("Fully independent data, check size of test.", {
  set.seed(730)
  num_sims <- 30
  n <- 20

  pval <- function(t) {
    X <- array(rnorm(n), dim=c(n))
    Y <- rnorm(n, 0, 1)
    result <- cdcv.univ.test(X,Y, M = 3, unbiased = TRUE, type = "circular")
    return(result$pCDCV)
  }
  pvals <- sapply(1:num_sims, pval)
  suppressWarnings({
    result <- wilcox.test(pvals, mu = 0.5)
    expect_gt(result$p.value, 0.05)
  })
})

test_that("Highly dependent data.", {
  set.seed(730)
  num_sims <- 30
  n <- 20

  Phis <- list(matrix(c(0,0.5,-0.5,0), 2, 2))
  Thetas <- list()
  Sigma <- matrix(c(1,0,0,1),2,2)

  pval <- function(t) {

    Z <- varma(n, Phis, Thetas, Sigma)
    X <- Z[,"X"]
    Y <- Z[,"Y"]
    result <- cdcv.univ.test(X,Y)
    return(result$pCDCV)
  }

  pvals <- sapply(1:num_sims, pval)
  suppressWarnings({
    result <- wilcox.test(pvals, alternative = "less", mu = 0.5)
    expect_lt(result$p.value, 0.05)
  })
})

