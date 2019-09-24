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
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
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
#' Simulates \eqn{n} points from \eqn{Step(X, Y) \in \mathbf{R}^d\times \mathbf{R}}{Step-Function(X, Y)} where:
#' \deqn{X \sim {U}\left(a, b\right)^d}{X ~ U(a, b)^d},
#' \deqn{Y = \mathbf{I}\left\{w^TX > 0\right\} + \kappa \epsilon N(0, 1)}{Y = I{w^TX > 0} + K*eps*N(0, 1)}
#' and \eqn{\kappa = 1\textrm{ if }d = 1, \textrm{ and 0 otherwise}}{K = 1 if d=1, and 0 otherwise} controls the noise for higher dimensions.
#'
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
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
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
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
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
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
#' Simumlates \eqn{n} points from \eqn{Spiral(X, Y) \in \mathbf{R}^d \times \mathbf{R}}{Spiral(X, Y)} where:
#' \eqn{X_i = U\, \textrm{cos}(\pi\, U)^d}{Xi = U*cos(pi*U)^d} if \code{i = d}, and \eqn{U\, \textrm{sin}(\pi U)\textrm{cos}^i(\pi U)}{sin(pi*U)*cos(pi*U)^i} otherwise
#' \deqn{Y = U\, \textrm{sin}(\pi\, U) + \epsilon p N(0, 1)}{Y = U*sin(pi*U) + eps*p*N(0, 1)}
#'
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
#'
#' @examples
#' library(mgc)
#' result  <- mgc.sims.spiral(n=100, d=10)  # simulate 100 samples in 10 dimensions
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
#' Simumlates \eqn{n} points from \eqn{Wshape(X, Y) \in \mathbf{R}^d \times \mathbf{R}}{Wshape(X, Y)} where:
#' \deqn{U \sim Bern(p)}{U ~ B(p)}
#' \deqn{X \sim Bern\left(p\right)^d + \epsilon N(0, I_d)}{X ~ B(p)^d + eps*N(0, I_d)}
#' \deqn{Y = (2U - 1)w^TX + \epsilon N(0, 1)}{Y = (2*U-1)w^T*X + 0.5*eps*N(0, 1)}
#'
#' For more details see the help vignette:
#' \code{vignette("sims", package = "mgc")}
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

#' Discriminability Linear Simulation
#'
#' A function to simulate multi-class data with a linear class-mean trend. The signal dimension is the dimension carrying all of the
#' between-class difference, and the non-signal dimensions are noise.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param rotate whether to apply a random rotation. Defaults to \code{TRUE}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
discr.sims.linear <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale
  S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*i*signal.lshift
  }), along=2)

  if (rotate) {
    res <- mgc.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- mgc.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=Y, mus=mus, Sigmas=S, priors=priors, simtype="Linear",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift, non.scale=non.scale,
                          rotate=rotate, class.equal=class.equal, ind=ind)))
}

#' Discriminability Exponential Simulation
#'
#' A function to simulate multi-class data with an Exponential class-mean trend.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions. The first dimension will be the signal dimension; the remainders noise.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param signal.lshift the location shift for the signal dimension between the classes. Defaults to \code{1}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param rotate whether to apply a random rotation. Defaults to \code{TRUE}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @author Eric Bridgeford
#' @export
discr.sims.exp <- function(n, d, K, signal.scale=1, signal.lshift=1, non.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S[1, 1] <- signal.scale; S[-c(1), -c(1)] <- non.scale
  S <- abind(lapply(1:K, function(i) {
    S
  }), along=3)
  mu <- c(1, rep(0, d-1))
  mu[1] <- signal.lshift
  mus <- abind(lapply(1:K, function(i) {
    mu*exp(i)*signal.lshift
  }), along=2)

  if (rotate) {
    res <- mgc.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
    rotate=res$Q
  }

  sim <- mgc.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)
  return(list(X=X, Y=Y, mus=mus, Sigmas=S, priors=priors, simtype="Linear",
              params=list(signal.scale=signal.scale, signal.lshift=signal.lshift,
                          non.scale=non.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}

#' Discriminability Spread Simulation
#'
#' A function to simulate data with the same mean that spreads as class id increases.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{1}.
#' @param rotate whether to apply a random rotation. Defaults to \code{TRUE}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(mgc)
#' sim <- discr.sims.fat_tails(100, 3, 2)
#' @author Eric Bridgeford
#' @export
discr.sims.fat_tails <- function(n, d, K, signal.scale=1, rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- signal.scale*diag(d)
  S <- abind(lapply(1:K, function(i) {
    i^2*S
  }), along=3)
  mus <- array(0, dim=c(d, K))

  if (rotate) {
    res <- mgc.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  sim <- mgc.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)

  if (ind) {
    X <- mgc.sims.sim_gmm(mus, S, n, priors=priors)$X
  }

  return(list(X=X, Y=Y, mus=mus, Sigmas=S, priors=priors, simtype="Fat Tails",
              params=list(signal.scale=signal.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}

#' Discriminability Cross Simulation
#'
#' A function to simulate data with the same mean that spreads as class id increases.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset.
#' @param signal.scale the scaling for the signal dimension. Defaults to \code{10}.
#' @param non.scale the scaling for the non-signal dimensions. Defaults to \code{1}.
#' @param mean.scale whether the magnitude of the difference in the means between the two classes.
#' If a mean scale is requested, \code{d} should be at least > \code{K}.
#' @param rotate whether to apply a random rotation. Defaults to \code{TRUE}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(mgc)
#' sim <- discr.sims.cross(100, 3, 2)
#' @author Eric Bridgeford
#' @export
discr.sims.cross <- function(n, d, K, signal.scale=10, non.scale=1, mean.scale=0,
                             rotate=FALSE, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)
  S <- diag(d)
  S <- abind(lapply(1:K, function(i) {
    S.class <- S
    S.class[i,i] <- signal.scale
    S.class[-c(i), -c(i)] <- non.scale
    S.class
  }), along=3)
  mus <- array(0, dim=c(d, K))
  if (mean.scale > 0) {
    if (d < K) {
      stop("If a mean scale is requested, d should be >= K.")
    }
    diag(mus) <- mean.scale
  }

  if (rotate) {
    res <- mgc.sims.random_rotate(mus, S)
    mus <- res$mus
    S <- res$S
  }

  sim <- mgc.sims.sim_gmm(mus, S, n, priors=priors)
  X <- sim$X; Y <- factor(sim$Y)

  if (ind) {
    X <- mgc.sims.sim_gmm(mus, S, n, priors=priors)$X
  }

  return(list(X=X, Y=Y, mus=mus, Sigmas=S, priors=priors, simtype="Logarithmic",
              params=list(signal.scale=signal.scale, rotate=rotate, class.equal=class.equal,
                          ind=ind)))
}
#' Discriminability Radial Simulation
#'
#' A function to simulate data with the same mean with radial symmetry as class id increases.
#'
#' @import abind
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param K the number of classes in the dataset.
#' @param er.scale the scaling for the error of the samples. Defaults to \code{0.1}.
#' @param r the radial spacing between each class. Defaults to \code{1}.
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @param ind whether to sample x and y independently. Defaults to \code{FALSE}.
#' @examples
#' library(mgc)
#' sim <- discr.sims.radial(100, 3, 2)
#' @author Eric Bridgeford
#' @export
discr.sims.radial <- function(n, d, K, er.scale=0.1, r=1, class.equal=TRUE, ind=FALSE) {
  priors <- gen.sample.labels(K, class.equal=class.equal)

  class <- 1:K
  # sample n points from 1:K for assignment of class label
  Y <- sample(class, size=n, replace=TRUE, prob=priors)
  # initialize X
  X <- array(NaN, dim=c(n, d))

  # loop over class labels that we are sampling from
  for (i in class[class %in% Y]) {
    if (i == 1) {
      f <- mgc.sims.2ball
    } else {
      f <- mgc.sims.2sphere
    }
    repvec <- Y == i
    pts <- do.call(f, list(sum(repvec), d, r=r*i, cov.scale=er.scale))
    X[repvec,] <- pts
  }

  if (ind) {
    Y <- sample(class, size=n, replace=TRUE, prob=priors)
  }
  return(list(X=X, Y=factor(Y), priors=priors, simtype="Radial",
              params=list(er.scale=er.scale, r=r, class.equal=class.equal,
                          ind=ind)))
}

#' A helper function for simulating sample labels
#' @param K the number of classes
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @keywords internal
gen.sample.labels <- function(K, class.equal=TRUE) {
  if (isTRUE(class.equal)) {
    priors <- array(1/K, dim=c(K)) # prior is just 1/K for all k
  } else {
    priors <- (1:K)/sum(1:K) # prior is k/sum(K) for all k
  }
  # return class priors
  return(priors)
}


#'
#' A helper function to generate a d-dimensional linear transformation matrix.
#' @param d the number of dimensions.
#' @return A \code{[d]} the coefficient vector.
#' @author Eric Bridgeford
#' @keywords internal
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
#' @keywords internal
gen.x.unif <- function(n, d, a=-1, b=1) {
  x <- array(runif(n=(n*d), min=a, max=b), dim=c(n, d))
  return(x)
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus \code{[d, K]} the mus for each class.
#' @param Sigmas \code{[d,d,K]} the Sigmas for each class.
#' @param n the number of examples.
#' @param priors \code{K} the priors for each class.
#' @return A list with the following:
#' \item{X}{\code{[n, d]} the simulated data.}
#' \item{Y}{\code{[n]} the labels for each data point.}
#' \item{priors}{\code{[K]} the priors for each class.}
#' @author Eric Bridgeford
#' @importFrom MASS mvrnorm
#' @keywords internal
mgc.sims.sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  labs <- sample(1:K, size=n, prob=priors, replace=TRUE)
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs, priors=priors))
}

#' Sample Random Rotation
#'
#' A helper function for estimating a random rotation matrix.
#' @importFrom stats rnorm
#' @param d dimensions to generate a rotation matrix for.
#' @return the rotation matrix
#' @author Eric Bridgeford
#' @keywords internal
mgc.sims.rotation <- function(d) {
  Q <- qr.Q(qr(array(rnorm(d*d), dim=c(d, d))))
  if (det(Q) < -.99) {
    Q[,1] <- -Q[,1]
  }
  return(Q)
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @param mus means per class.
#' @param Sigmas covariances per class.
#' @param Q rotation to use, if any
#' @author Eric Bridgeford
#' @keywords internal
mgc.sims.random_rotate <- function(mus, Sigmas, Q=NULL) {
  dimm <- dim(mus)
  K <- dimm[2]
  d <- dim(mus)[1]
  if (is.null(Q)) {
    Q <- mgc.sims.rotation(d)
  } else if (!isTRUE(all.equal(dim(Q), c(d, d)))) {
    stop(sprintf("You have specified a rotation matrix with dimensions (%d, %d), but should be (%d, %d).", dim(Q)[1], dim(Q)[2], d, d))
  }

  for (i in 1:K) {
    mus[,i] <- Q %*% mus[,i,drop=FALSE]
    Sigmas[,,i] <- Q %*% Sigmas[,,i] %*% t(Q)
  }
  return(list(mus=mus, S=Sigmas, Q=Q))
}

#' Sample from Unit 2-Ball
#'
#' Sample from the 2-ball in d-dimensions.
#'
#' @importFrom MASS mvrnorm
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{NaN},
#' which has no noise.
#' @return the points sampled from the ball, as a \code{[n, d]} array.
#' @examples
#' library(mgc)
#' # sample 100 points from 3-d 2-ball with radius 2
#' X <- mgc.sims.2ball(100, 3, 2)
#' @author Eric Bridgeford
#' @export
mgc.sims.2ball <- function(n, d, r=1, cov.scale=0) {
  Y <- matrix(mvrnorm(n=n, mu=array(0, dim=c(d, 1)), Sigma=diag(d)), nrow=n)
  u <- runif(n)
  r <- r * u^(1/d)
  X <- r * Y/sqrt(apply(Y^2, 1, sum)) + matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d)), nrow=n)
  return(X)
}

#' Sample from Unit 2-Sphere
#'
#' Sample from the 2-sphere in d-dimensions.
#'
#' @importFrom MASS mvrnorm ginv
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{0},
#' which has no noise.
#' @return the points sampled from the sphere, as a \code{[n, d]} array.
#' @examples
#' library(mgc)
#' # sample 100 points from 3-d 2-sphere with radius 2
#' X <- mgc.sims.2sphere(100, 3, 2)
#' @author Eric Bridgeford
#' @export
mgc.sims.2sphere <- function(n, d, r, cov.scale=0) {
  u <- matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=diag(d)), nrow=n)
  unorm <- diag(sqrt(apply(u^2, 1, sum)))
  pts <- r*(ginv(unorm) %*% u)
  pts <- pts + matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d)), nrow=n)
  return(pts)
}
