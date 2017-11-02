#' A simple demo function to illustrate MGC on test data.
#'
#' @export
mgc.run_demo <- function(x, y){
  set.seed(12345)
  mgc.sample(x=x, y=y)
}
#' A simple demo function to illustrate Discriminability on 2 realizations of 2 random variables.
#'
#' @export
discr.run_demo <- function() {
  n <- 2
  d <- 5
  s <- 2
  dat <- array(NaN, dim=c(n*s, d))
  mu1 <- 0
  mu2 <- 2
  sd1 <- .2
  sd2 <- .2
  dat[1,] <- rnorm(n=d, mean=mu1, sd=sd1)
  dat[2,] <- rnorm(n=d, mean=mu1, sd=sd1)
  dat[3,] <- rnorm(n=d, mean=mu2, sd=sd2)
  dat[4,] <- rnorm(n=d, mean=mu2, sd=sd2)

  print(sprintf("Parameters of first random variable: mu=%.2f, sd=%.2f", mu1, sd1))
  print(sprintf("Parameters of second random variable: mu=%.2f, sd=%.2f", mu2, sd2))
  D <- discr.distance(dat)
  print(sprintf("Discriminability statistic for 2 observations per random variable: %.2f", discr.discr(D, c(1,1,2,2))))
}
