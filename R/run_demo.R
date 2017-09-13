#' A simple demo function to illustrate MGC on linear and quadratic relationships at size 20.
#'
#' @export
mgc.run_demo <- function(){
  n=20;
  a=rnorm(n, mean = 0, sd = 1)
  b=3*a; # linear relationship
  rep=1000;
  mgc1=mgc.test(a,b,rep) # mgc test with p-value 0 and optimal scale always at the maximal (400)
  cat('Under Linear Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are: ', mgc1$statMGC, ',', mgc1$pMGC,',',mgc1$optimalScale,'\n')

  b=a^2; # quadratic relationship
  mgc2=mgc.test(a,b,rep) #mgc test with p-value almost 0, and local optimal scale
  cat('Under Quadratic Dependency of 20 observations, MGC Statistic, p-value, and optimal scales (matrix single index) are: ', floor(mgc2$statMGC*100)/100, ',', mgc2$pMGC,',',mgc2$optimalScale,'\n')
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
