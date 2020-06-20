context("discr")
require(MASS)
require(abind)

sim_gmm <- function(mus, Sigmas, n) {
  K <- dim(mus)[2]
  ni <- round(n/K)
  labs <- c(sapply(1:K, function(k) rep(k, ni)))
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]),
                USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs))
}

## Linear Signal Difference
# a simulation where classes are linearly distinguishable
# 2 classes
sim.linear_sig <- function(n, d, sigma=2) {
  S <- diag(d)
  S[1, 1] <- 1
  S[-c(1), -c(1)] <- 1
  S2 <- S*sigma  # sample 2 has different covariance in signal dimension
  mus=cbind(rep(0, d), c(1, rep(0, d-1))) # with the same mean signal shift between the classes
  # sample 1 should be more discriminable than sample 2
  samp1 <- sim_gmm(mus=mus, Sigmas=abind(S, S, along=3), n)
  samp2 <- sim_gmm(mus=mus, Sigmas=abind(S2, S2, along=3), n)
  return(list(X1=samp1$X, X2=samp2$X, Y=samp1$Y))
}

test_that("RDF - 2 Class, d==1", {
  X <- as.matrix(c(0, 1, 2, 3)); Y <- c(1, 1, 2, 2)
  D <- mgc.distance(X)
  drdf <- discr.rdf(D, Y)
  expect_equal(as.numeric(drdf), c(1, .75, .75, 1))
})

test_that("Discriminability - 2 Class, d==1", {
  X <- as.matrix(c(0, 1, 2, 3)); Y <- c(1, 1, 2, 2)
  dstat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)
  expect_equal(dstat$discr, mean(c(1, .75, .75, 1)))

  X <- as.matrix(c(0, 0, 2, 2)); Y <- c(1, 1, 2, 2)
  dstat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)
  expect_equal(dstat$discr, mean(c(1, 1, 1, 1)))
  expect_equal(as.numeric(dstat$rdf), c(1, 1, 1, 1))
})

test_that("Discriminability - 2 Class, d > 1; large n", {
  n=50; d=3
  # large class difference
  sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=20); X <- sim$X; Y <- sim$Y
  dstat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)
  expect_equal(dstat$discr, 1)

  # small class difference
  sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=.05); X <- sim$X; Y <- sim$Y
  dstat.dat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)

  # large class difference is more discriminable
  expect_lte(dstat.dat$discr, dstat$discr)
})

test_that("Discriminability - using data and distance function manually produces same result", {
  n=100; d=3

  # small class difference with distance matrix vs data produces same result
  sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=2); X <- sim$X; Y <- sim$Y
  dstat.dat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)

  D <- mgc.distance(X)
  dstat.dist <- discr.stat(D, Y, is.dist = TRUE)
  expect_equal(dstat.dat$discr, dstat.dist$discr)
})

test_that("One Sample Test is Valid", {
  n = 100; d=2; nsim=10; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(lapply(1:nsim, function(i) {
    # no true class difference
    set.seed(seed.idx[i])
    sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=0); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(discr.test.one_sample(X, Y, nperm=50)$p.value < alpha)
  }), use.names=FALSE)
  # check power is near alpha
  expect_lte(abs(mean(res) - alpha), 0.1)
})

test_that("One Sample Test Detects Relationship", {
  n = 50; d=3; nsim=5; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(lapply(1:nsim, function(i) {
    # substantial true class difference
    set.seed(seed.idx[i])
    sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=4); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(discr.test.one_sample(X, Y, nperm=20)$p.value < alpha)
  }), use.names=FALSE)
  # check power is near 1
  expect_lte(abs(mean(res) - 1), 0.1)
})

test_that("Two Sample Test is Valid and Detects Relationship", {
  n = 50; d=2; nsim=10; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  sim.opts <- list(list(alt="greater", sigma=20, outcome=1),
                   list(alt="neq", sigma=30, outcome=1),
                   list(alt="less", sigma=0.01, outcome=1))
  # test all cases of alternatives that can be specified
  sapply(sim.opts, function(sim.opt) {
    res <- unlist(lapply(1:nsim, function(i) {
      # true class difference, where dataset 1 more discriminable than dataset 2
      set.seed(seed.idx[i])
      sim <- sim.linear_sig(n, d, sigma=sim.opt$sigma)
      set.seed(seed.idx[i])
      pval <- discr.test.two_sample(sim$X1, sim$X2, sim$Y, alt=sim.opt$alt, nperm=20)$p.value
      return(pval < alpha)
    }), use.names=FALSE)
    # check power accordingly
    expect_lte(abs(mean(res) - sim.opt$outcome), .2)
  })
})
