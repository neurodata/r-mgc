#' Linear Simulation
#'
#' A function for Generating a linear simulation.
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
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbb{R}^d \cross \mathbb{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
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
  xs <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- xs%*%w + kappa*eps*nu  # y = xA + nu
  if (ind) {
    xs <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Exponential Simulation
#'
#' A function for Generating an exponential simulation.
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
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbb{R}^d \cross \mathbb{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
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
  xs <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- exp(xs%*%w) + kappa*eps*nu  # y = exp(xA) + nu
  if (ind) {
    xs <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Cubic Simulation
#'
#' A function for Generating a cubic simulation.
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
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{w_i = \frac{1}{i}}{w[i] = 1/i} is a weight-vector that scales with the dimensionality.
#' Simulates \eqn{n} points from \eqn{Linear(X, Y) \in  \mathbb{R}^d \cross \mathbb{R}}{Linear(X, Y)}, where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d}
#' \deqn{Y = c_3\left(w^TX - s\right)^3 + c_2\left(w^TX - s\right)^2 + c_1\left(w^TX - s\right) + \kappa \epsilon}{Y = c[3](w^TX - s)^3 + c[2](w^TX - s)^2 + c[1](w^TX - s) + \kappa \epsilon}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.cubic(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.cubic <- function(n, d, eps=80, ind=FALSE, a=-1, b=1, c=c(-12, 48, 128), s=1/3) {
  xs <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  xw = xs%*%w
  ys <- c[3]*(xw - s)^3 + c[2]*(xw - s)^2 + c[1]*(xw - s) + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' Joint Normal Simulation
#'
#' A function for Generating a joint-normal simulation.
#'
#' @importFrom MASS mvrnorm
#' @param n the number of samples for the simulation.
#' @param d the number of dimensions for the simulation setting.
#' @param eps the noise level for the simulation. Defaults to \code{0.5}.
#' @return a list containing the following:
#' \item{\code{X}}{\code{[n, d]} the data matrix with \code{n} samples in \code{d} dimensions.}
#' \item{\code{Y}}{\code{[n]} the response array.}
#'
#' @section Details:
#' Given: \eqn{\rho = \frac{1}{2}d}{r = 1/2*d}, \eqn{I_d}{Id} is the identity matrix of size \eqn{d \times d}{dxd}, \eqn{J_d}{Jd} is the matrix of ones of size \eqn{d \times d}{dxd}.
#' Simulates \eqn{n} points from \eqn{Joint-Normal(X, Y) \in  \mathbb{R}^d \cross \mathbb{R}^d}{Joint-Normal(X, Y)}, where:
#' \deqn{(X, Y) \sim \mathcal{N}(0, \Sigma)}{(X, Y) ~ N(0, E)},
#' \deqn{\Sigma = \begin{bmatrix}I_d & \rho J_d \\ \rho J_d & (1 + \epsilon\kappa)I_d\end{bmatrix}}{E = [Id, r*Jd; r*Jd, (1+eps*K)*Id]}
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
  y = samp[, (d+1):(2*d)] + kappa*eps*nu
  x = samp[, 1:d]
  return(list(X=x, Y=y))
}

#' Step Function Simulation
#'
#' A function for Generating a step function simulation.
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
#' Simulates \eqn{n} points from \eqn{Step-Function(X, Y) \in \mathbb{R}^d\times \mathbb{R}}{Step-Function(X, Y)} where:
#' \deqn{X \sim \mathcal{U}\left(a, b\right)^d}{X ~ U(a, b)^d},
#' \deqn{Y = \mathbb{I}\left\{w^TX > 0\right\} + \epsilon}{Y = I{w^TX > 0} + eps}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.step(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.step <- function(n, d, eps=1, ind=FALSE, a=-1, b=1) {
  xs <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- as.numeric(xs%*%w > 0) + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' Quadratic Simulation
#'
#' A function for Generating a quadratic simulation.
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
#' Simulates \code{n} points from \eqn{Quadratic(X, Y) \in \mathbb{R}^d \times \mathbb{R}}{Quadratic(X, Y)} where:
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d},
#' \deqn{Y = (w^TX)^2 + \kappa\epsilon}{Y = (w^TX)^2 + K*eps}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.quad(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.quad <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- (xs%*%w)^2 + eps*kappa*nu
  if (ind) {
    xs <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' W Shaped Simulation
#'
#' A function for Generating a W-shaped simulation.
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
#' Simumlates \eqn{n} points from \eqn{W-shape(X, Y) \in \mathbb{R}^d \times \mathbb{R}}{W-shape(X, Y)} where:
#' \deqn{U \sim \mathcal{U}(a, b)^d}{U ~ U(a, b)^d},
#' \deqn{X \sim \mathcal{U}(a, b)^d}{X ~ U(a, b)^d},
#' \deqn{Y = \left[\left((w^TX)^2 - \frac{1}{2}\right)^2 + \frac{w^TU}{500}\right] + \kappa \epsilon}{Y = [((w^TX)^2 - 1/2)^2 + w^TU/500] + K*eps}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.wshape(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.wshape <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  x <- gen.x(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- 4*(((x%*%w)^2 - 1/2)^2 + u%*%w/500) + kappa*eps*nu
  if (ind) {
    x <- gen.x(n, d, a=a, b=b)
  }
  return(list(X=x, Y=y))
}

#' Spiral Simulation
#'
#' A function for Generating a spiral simulation.
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
#' Given: \eqn{U \sim \mathcal{U}(a, b)}{U ~ U(a, b)} a random variable.
#' \deqn{X_i = \begin{cases}U\textrm{sin}(\pi U)cos^d(\pi U) & i < d \\ \end{cases}}{Xi = U*cos(pi*U)^d if i = d, and Usin(pi*U)*cos(pi*U)^d otherwise}
#' \deqn{Y = U\textrm{sin}(\pi U) + \epsilon p}{Y = U*sin(pi*U) + eps*p}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.wshape(n=100, d=10)  # simulate 100 samples in 10 dimensions
#' X <- result$X; Y <- result$Y
#' @author Eric Bridgeford
#' @export
mgc.sims.spiral <- function(n, d, eps=0.4, a=0, b=5) {
  u <- gen.x(n, 1, a=a, b=b)
  x <- array(cos(pi*u), dim=c(n, d))
  y <- u*sin(pi*u)
  for (i in 1:(d-1)) {
    x[, i] <- y*ru[, i]^i
  }
  x[, d] <- u*ru[, d]
  kappa <- as.numeric(d == 1)
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

#' Uncorrelated Bernoulli Simulation
#'
#' A function for Generating an uncorrelated bernoulli simulation.
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
#' Simulates \eqn{n} points from \eqn{UBern(X, Y) \in  \mathbb{R}^d \cross \mathbb{R}}{UBern(X, Y)}, where:
#' \deqn{U \sim \mathcal{B}(p)}{U ~ B(p)}
#' \deqn{X \sim \mathcal{B}\left(p\right)^d + \epsilon\mathcal{(0, I_d)}}{X ~ B(p)^d + eps*N(0, I_d)}
#' \deqn{Y = (2U - 1)w^TX + \epsilon}{Y = (2*U-1)w^T*X + 0.5*eps}
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
  #nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  w <- gen.coefs(d)
  Y <- array(NaN, dim=c(n, 1))
  nu_e2 <- rnorm(n, mean=0, sd=1)
  for (i in 1:n) {
    Y[i] <- (2*U[i] - 1)*w^T %*% X[i,] + eps*nu_e2[i]
  }
  return(list(X=X, Y=Y))
}

#' A helper function to generate a d-dimensional linear transformation matrix.
gen.coefs <- function(d) {
  A = as.array(1/1:d, dim=c(d, 1))
  return(A)
}

#' A helper function to generate n samples of a d-dimensional uniform vector.
gen.x <- function(n, d, a=-1, b=1) {
  x <- as.array(replicate(d, runif(n, a, b)))
  return(t(t(x)))
}
