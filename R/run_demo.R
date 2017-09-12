#' A simple demo function to illustrate MGC on linear and quadratic relationships at size 20.
#'
#' @export
run_demo <- function(){
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
