#' Linear Simulation
#'
#' A function for Generating a linear simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{1}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{1}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbf{R}^d \times \mathbf{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim {U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = w^TX + \kappa \epsilon}{Y = w^T X + \kappa \epsilon}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.linear(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.linear <- function(n, d, eps=1, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- xs%*%w + kappa*eps*nu  # y = xA + nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Exponential Simulation
#'
#' A function for Generating an exponential simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{10}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{0}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{3}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbf{R}^d \times \mathbf{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim {U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = e^{w^TX} + \kappa \epsilon}{Y = exp(w^T X) + \kappa \epsilon}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.exp(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.exp <- function(n, d, eps=10, ind=FALSE, a=0, b=3) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- exp(xs%*%w) + kappa*eps*nu  # y = exp(xA) + nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Cubic Simulation
#'
#' A function for Generating a cubic simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{80}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the range of the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the range  of the data matrix. Defaults to \code{1}.
#' @param c.coef the coefficients for the cubic function, where the first value is the first order coefficient, the second value the quadratic coefficient, and the third the cubic coefficient. Defaults to \code{c(-12, 48, 128)}.
#' @param s the scaling for the center of the cubic. Defaults to \code{1/3}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbf{R}^d \times \mathbf{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim {U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = c_3\left(w^TX - s\right)^3 + c_2\left(w^TX - s\right)^2 + c_1\left(w^TX - s\right) + \kappa \epsilon}{Y = c[3](w^TX - s)^3 + c[2](w^TX - s)^2 + c[1](w^TX - s) + \kappa \epsilon}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.cubic(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.cubic <- function(n, d, eps=80, ind=FALSE, a=-1, b=1, c.coef=c(-12, 48, 128), s=1/3) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  xw = xs%*%w
  ys <- c.coef[3]*(xw - s)^3 + c.coef[2]*(xw - s)^2 + c.coef[1]*(xw - s) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Joint Normal Simulation
#'
#' A function for Generating a joint-normal simulation.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{\rho = \frac{1}{2}d}{r = 1/2*d}, \eqn{I_d}{Id} is the identity matrix of size \eqn{d \times d}{dxd}, \eqn{J_d}{Jd} is the matrix of ones of size \eqn{d \times d}{dxd}.
#' Simulates \eqn{n} points from \eqn{Joint-Normal(X, Y) \in  \mathbf{R}^d \times \mathbf{R}^d}{Joint-Normal(X, Y)}, where:
#' \deqn{(X, Y) \sim {N}(0, \Sigma)}{(X, Y) ~ N(0, E)},
#' \deqn{\Sigma = \left[I_d, \rho J_d; \rho J_d , (1 + \epsilon\kappa)I_d\right]}{E = [Id, r*Jd; r*Jd, (1+eps*K)*Id]}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.joint(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.joint <- function(n, d, eps=0.5) {
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  rho = 1/(2*d)
  Sig <- diag(2*d)
  Sig[(d+1):(2*d), (1:d)] = rho
  Sig[(1:d), (d+1):(2*d)] = rho
  samp = mvrnorm(n=n, mu=array(0, dim=c(2*d, 1)), Sigma=Sig)
  y = samp[, (d+1):(2*d), drop=FALSE] + kappa*eps*nu
  x = samp[, 1:d, drop=FALSE]
  return(list(X=x, Y=y))
}

#' Step Function Simulation
#'
#' A function for Generating a step function simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{1}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the data matrix. Defaults to \code{-1}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Step-Function(X, Y) \in \mathbf{R}^d\times \mathbf{R}}{Step-Function(X, Y)} where:
#' \deqn{X \sim {U}\left(a, b\right)^d}{X ~ U(a, b)^d},
#' \deqn{Y = \mathbf{I}\left\{w^TX > 0\right\} + \kappa \epsilon N(0, 1)}{Y = I{w^TX > 0} + K*eps*N(0, 1)}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.step(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.step <- function(n, d, eps=1, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- as.numeric(xs%*%w > 0) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' Quadratic Simulation
#'
#' A function for Generating a quadratic simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the data matrix. Defaults to \code{-1}.
#' @param b the upper limit for the data matrix. Defaults to \code{1}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \code{n} points from \eqn{Quadratic(X, Y) \in \mathbf{R}^d \times \mathbf{R}}{Quadratic(X, Y)} where:
#' \deqn{X \sim {U}(a, b)^d}{X ~ U(a, b)^d},
#' \deqn{Y = (w^TX)^2 + \kappa\epsilon N(0, 1)}{Y = (w^TX)^2 + K*eps*N(0, 1)}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.quad(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.quad <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- (xs%*%w)^2 + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' W Shaped Simulation
#'
#' A function for Generating a W-shaped simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @param a the lower limit for the data matrix. Defaults \code{-1}.
#' @param b the upper limit for the data matrix. Defaults to \code{1}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simumlates \eqn{n} points from \eqn{W-shape(X, Y) \in \mathbf{R}^d \times \mathbf{R}}{W-shape(X, Y)} where:
#' \deqn{U \sim {U}(a, b)^d}{U ~ U(a, b)^d},
#' \deqn{X \sim {U}(a, b)^d}{X ~ U(a, b)^d},
#' \deqn{Y = \left[\left((w^TX)^2 - \frac{1}{2}\right)^2 + \frac{w^TU}{500}\right] + \kappa \epsilon N(0, 1)}{Y = [((w^TX)^2 - 1/2)^2 + w^TU/500] + K*eps*N(0, 1)}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.wshape(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.wshape <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  x <- gen.x.unif(n, d, a=a, b=b)
  u <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- 4*(((x%*%w)^2 - 1/2)^2 + u%*%w/500) + kappa*eps*nu
  if (ind) {
    x <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=x, Y=y))
}

#' Spiral Simulation
#'
#' A function for Generating a spiral simulation.
#'
#' @importFrom stats rnorm
#'
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @param a the lower limit for the data matrix. Defaults \code{-1}.
#' @param b the upper limit for the data matrix. Defaults to \code{1}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{U \sim U(a, b)}{U ~ U(a, b)} a random variable.
#' \eqn{X_i = {U\textrm{cos}(\pi U)^d}{Xi = U*cos(pi*U)^d} if \code{i = d}, and eqn{U\textrm{sin}(\pi U)cos^i(\pi U)}{sin(pi*U)*cos(pi*U)^i} otherwise
#' \deqn{Y = U\textrm{sin}(\pi U) + \epsilon p N(0, 1)}{Y = U*sin(pi*U) + eps*p*N(0, 1)}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.wshape(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.spiral <- function(n, d, eps=0.4, a=0, b=5) {
  u <- gen.x.unif(n, 1, a=a, b=b)
  x <- array(cos(pi*u), dim=c(n, d))
  y <- u*sin(pi*u)
  if (d > 1) {
    for (i in 1:(d-1)) {
      x[, i] <- y*x[, i, drop=FALSE]^i
    }
  }
  x[, d] <- u*x[, d, drop=FALSE]
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

#' Uncorrelated Bernoulli Simulation
#'
#' A function for Generating an uncorrelated bernoulli simulation.
#'
#' @importFrom stats rnorm rbinom
#' @importFrom MASS mvrnorm
#'
#' @importFrom MASS mvrnorm
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @param p the bernoulli probability.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{UBern(X, Y) \in  \mathbf{R}^d \times \mathbf{R}{UBern(X, Y)}, where:
#' \deqn{U \sim Bern(p)}{U ~ B(p)}
#' \deqn{X \sim Bern\left(p\right)^d + \epsilon N(0, I_d)}{X ~ B(p)^d + eps*N(0, I_d)}
#' \deqn{Y = (2U - 1)w^TX + \epsilon N(0, 1)}{Y = (2*U-1)w^T*X + 0.5*eps*N(0, 1)}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.ubern(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.ubern <- function(n, d, eps=0.5, p=0.5) {
  U = rbinom(n=n, size=1, prob=p)
  nu_e1 <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))  # gaussian noise
  X = array(rbinom(n=d*n, size=1, prob=p), dim=c(n, d)) + eps*nu_e1
  w <- gen.coefs(d)
  Y <- array(NaN, dim=c(n, 1))
  nu_e2 <- rnorm(n, mean=0, sd=1)
  for (i in 1:n) {
    Y[i] <- (2*U[i] - 1)*w^T %*% X[i,] + eps*nu_e2[i]
  }
  return(list(X=X, Y=Y))
}

#' A helper function to generate a d-dimensional linear transformation matrix.
#' @param d the number of dimensions.
#' @return A \code{[d]} the coefficient vector.
#' @author Eric Bridgeford
gen.coefs <- function(d) {
  A = as.array(1/1:d, dim=c(d, 1))
  return(A)
}

#' A helper function to generate n samples of a d-dimensional uniform vector.
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param a the lower limit.
#' @param b the upper limit.
#' @param x \code{[n, d]} the simulated data matrix.
#' @author Eric Bridgeford
#' @importFrom stats runif
gen.x.unif <- function(n, d, a=-1, b=1) {
  x <- array(runif(n=(n*d), min=a, max=b), dim=c(n, d))
  return(x)
}
