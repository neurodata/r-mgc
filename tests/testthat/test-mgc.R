test_that("MGC Detects Global Scale when Global Scale True", {
  n = 100; d=5; nsim=50; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  result <- parallel::mclapply(1:nsim, function(i) {
    # true, linear dependence
    set.seed(seed.idx[i])
    sim <- mgc.sims.linear(n=n, d=d, eps=0.5, ind=FALSE); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    res <- mgc.test(X, Y, nperm=100)
    return(list(test=res$p.value < alpha, scale=res$optimalScale))
  }, mc.cores=parallel::detectCores() - 1)
  # check power is near alpha
  expect_lt(abs(mean(sapply(result, function(res) res$test)) - 0.9), 0.1)
  # check optimal scale is no more than a radius of 5 from c(n, n)
  expect_lt(norm(Reduce("+", lapply(result, function(res) as.numeric(res$scale)))/nsim - c(n, n), type="2"), 50)
})

test_that("MGC is Valid", {
  n = 100; d=5; nsim=50; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(parallel::mclapply(1:nsim, function(i) {
    # no true class difference
    set.seed(seed.idx[i])
    sim <- mgc.sims.linear(n=n, d=d, eps=1, ind=TRUE); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(mgc.test(X, Y, nperm=100)$p.value < alpha)
  }, mc.cores=parallel::detectCores() - 1), use.names=FALSE)
  # check power is near alpha
  expect_lt(abs(mean(res) - alpha), 0.1)
})

test_that("MGC Works with Distance and Non-Distance Matrices for X and Y", {
  n = 100; d=5; nsim=50; alpha=0.1
  sim <- mgc.sims.linear(n=n, d=d, eps=1, ind=FALSE); X <- sim$X; Y <- sim$Y
  DX <- mgc.distance(sim$X); DY <- mgc.distance(sim$Y)
  test.xy <- mgc.test(X, Y, nperm=100, no_cores=parallel::detectCores() - 1)
  test.dx <- mgc.test(DX, Y, nperm=100, is.dist.X=TRUE, no_cores = parallel::detectCores() - 1)
  test.dy <- mgc.test(X, DY, nperm=100, is.dist.Y=TRUE, no_cores = parallel::detectCores() - 1)
  test.dxy <- mgc.test(DX, DY, nperm=100, is.dist.X=TRUE, is.dist.Y=TRUE, no_cores = parallel::detectCores() - 1)

  expect_true(all.equal(test.xy$p.value, test.dx$p.value, test.dy$p.value, test.dxy$p.value))
  expect_true(all.equal(test.xy$stat, test.dx$stat, test.dy$stat, test.dxy$stat))
})
