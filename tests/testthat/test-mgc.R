test_that("MGC Detects Global Scale when Global Scale True", {
  n = 50; d=5; nsim=10; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  result <- lapply(1:nsim, function(i) {
    # true, linear dependence
    set.seed(seed.idx[i])
    sim <- mgc.sims.linear(n=n, d=d, eps=0.5, ind=FALSE); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    res <- suppressWarnings(mgc.test(X, Y, nperm=20))
    return(list(test=res$p.value < alpha, scale=res$optimalScale))
  })
  # check power is near alpha
  expect_lte(abs(mean(sapply(result, function(res) res$test)) - 0.9), 0.1)
  # check optimal scale is no more than a radius of 5 from c(n, n)
  expect_lte(norm(Reduce("+", lapply(result, function(res) as.numeric(res$scale)))/nsim - c(n, n), type="2"), 50)
})

test_that("MGC is Valid", {
  n = 50; d=2; nsim=15; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(lapply(1:nsim, function(i) {
    # no true class difference
    set.seed(seed.idx[i])
    sim <- mgc.sims.linear(n=n, d=d, eps=.5, ind=TRUE); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(suppressWarnings(mgc.test(X, Y, nperm=10)$p.value < alpha))
  }), use.names=FALSE)
  # check power is near alpha
  expect_lte(abs(mean(res) - alpha), 0.15)
})

test_that("MGC Works with Distance and Non-Distance Matrices for X and Y", {
  n = 50; d=5; alpha=0.1
  sim <- mgc.sims.linear(n=n, d=d, eps=.5, ind=FALSE); X <- sim$X; Y <- sim$Y
  DX <- mgc.distance(sim$X); DY <- mgc.distance(sim$Y)
  test.xy <- suppressWarnings(mgc.test(X, Y, nperm=10, no_cores=1))
  test.dx <- suppressWarnings(mgc.test(DX, Y, nperm=10, is.dist.X=TRUE, no_cores = 1))
  test.dy <- suppressWarnings(mgc.test(X, DY, nperm=10, is.dist.Y=TRUE, no_cores = 1))
  test.dxy <- suppressWarnings(mgc.test(DX, DY, nperm=10, is.dist.X=TRUE, is.dist.Y=TRUE, no_cores = 1))

  expect_true(all.equal(test.xy$p.value, test.dx$p.value, test.dy$p.value, test.dxy$p.value, tolerance = 0.05))
  expect_true(all.equal(test.xy$stat, test.dx$stat, test.dy$stat, test.dxy$stat))
})
