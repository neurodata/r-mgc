#' Linear Simulation
#'
#' A function for Generating a Linear Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{1}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{1}.
#' @return a list containing the following:
#' \itemize{
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#' }
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathcal{R}^d \cross \mathcal{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = w^TX + \kappa \epsilon}{Y = w^T X + \kappa \epsilon}
#' and \eqn{\kappa = 1 if d = 1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.linear(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.linear <- function(n, d, eps=1, ind=FALSE, a=-1, b=-1) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- xs%*%w + kappa*eps*nu  # y = xA + nu
  if (ind) {
    xs <- gen.x(n, d)
  }
  return(list(X=xs, Y=ys))
}

#' Exponential Simulation
#'
#' A function for Generating an Exponential Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{10}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{0}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{3}.
#' @return a list containing the following:
#' \itemize{
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#' }
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathcal{R}^d \cross \mathcal{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = e^{w^TX} + \kappa \epsilon}{Y = exp(w^T X) + \kappa \epsilon}
#' and \eqn{\kappa = 1 if d = 1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.exp(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.exp <- function(n, d, eps=10, ind=FALSE, a=0, b=3) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- exp(xs%*%w) + kappa*eps*nu  # y = exp(xA) + nu
  if (ind) {
    xs <- gen.x(n, d)
  }
  return(list(X=xs, Y=ys))
}

#' Cubic Simulation
#'
#' A function for Generating a Cubic Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{80}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{1}.
#' @param c the coefficients for the cubic function, where the first value is the first order coefficient, the second value the quadratic coefficient, and the third the cubic coefficient. Defaults to \code{c(-12, 48, 128)}.
#' @param s the scaling for the center of the cubic. Defaults to \code{1/3}.
#' @return a list containing the following:
#' \itemize{
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#' }
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathcal{R}^d \cross \mathcal{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = c_3\left(w^TX - s\right)^3 + c_2\left(w^TX - s\right)^2 + c_1\left(w^TX - s\right) + \kappa \epsilon}{Y = c[3](w^TX - s)^3 + c[2](w^TX - s)^2 + c[1](w^TX - s) + \kappa \epsilon}
#' and \eqn{\kappa = 1 if d = 1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.cubic(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.cubic <- function(n, d, eps=80, ind=FALSE, a=-1, b=-1, c=c(-12, 48, 128), s=1/3) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  xw = xs%*%w
  ys <- c[3]*(xw - s)^3 + c[2]*(xw - s)^2 + c[1]*(xw - s) + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d)
  }
  return(list(X=xs, Y=ys))
}

#' Joint Normal
#'
#' A function for Generating a Joint Normal Simulation.
#'
#' @import MASS
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps=0.5 the noise level for the simulation.
#' @return X [n, d] the data matrix.
#' @return Y [n] the response array.
#' @author Eric Bridgeford
#' @export
mgc.sims.joint <- function(n, d, eps=0.5) {
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  rho = 1/(2*d)
  Sig <- diag(2*d)
  Sig[(d+1):(2*d), (1:d)] = rho
  Sig[(1:d), (d+1):(2*d)] = rho
  samp = MASS::mvrnorm(n=n, mu=array(0, dim=c(2*d, 1)), Sigma=Sig)
  y = samp[, (d+1):(2*d)] + kappa*eps*nu
  x = samp[, 1:d]
  return(list(X=x, Y=y))
}

#' Step Function
#'
#' A function for Generating a Step Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps=1 the noise level for the simulation.
#' @param ind=FALSE whether to sample x and y independently.
#' @param a=-1 the lower limit for the data matrix.
#' @param b=1 the upper limit for the data matrix.
#' @return X [n, d] the data matrix.
#' @return Y [n] the response array.
#' @author Eric Bridgeford
#' @export
mgc.sims.step <- function(n, d, eps=1, ind=FALSE, a=-1, b=-1) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- as.numeric(x%*%w > 0) + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d)
  }
  return(list(X=x, Y=y))
}

#' Quadratic Function
#'
#' A function for Generating a Quadratic Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps=0.5 the noise level for the simulation.
#' @param ind=FALSE whether to sample x and y independently.
#' @param a=-1 the lower limit for the data matrix.
#' @param b=1 the upper limit for the data matrix.
#' @return X [n, d] the data matrix.
#' @return Y [n] the response array.
#' @author Eric Bridgeford
#' @export
mgc.sims.quad <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=-1) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- (x%*%w)^2 + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d)
  }
  return(list(X=x, Y=y))
}

#' W Shaped Function
#'
#' A function for Generating a W Simulation.
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps=0.5 the noise level for the simulation.
#' @param ind=FALSE whether to sample x and y independently.
#' @param a=-1 the lower limit for the data matrix.
#' @param b=1 the upper limit for the data matrix.
#' @return X [n, d] the data matrix.
#' @return Y [n] the response array.
#' @author Eric Bridgeford
#' @export
mgc.sims.wshape <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=-1) {
  x <- gen.x(n, d, a=a, b=b)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  nu <- rnorm(dim(x)[1], mean=0, sd=eps)  # gaussian noise
  y <- 4*(((x%*%w)^2 - 1/2)^2 + u%*%w/500)+ kappa*eps*nu
  return(list(X=x, Y=y))
}

#' A helper function to generate a d-dimensional linear transformation matrix.
gen.coefs <- function(d) {
  A = as.array(1/1:d, dim=c(d, 1))
  return(A)
}

#' A helper function to generate n samples of a d-dimensional uniform vector.
gen.x <- function(n, d, a=-1, b=-1) {
  x <- as.array(replicate(d, runif(n, a, b)))
  return(t(t(x)))
}
