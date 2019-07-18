context("discr")

test_that("RDF - 2 Class, d==1", {
  X <- as.matrix(c(0, 1, 2, 3)); Y <- c(1, 1, 2, 2)
  D <- discr.distance(X)
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

test_that("Discriminability - 2 Class, d > 1", {
  n=100; d=10
  # large class difference
  sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=10); X <- sim$X; Y <- sim$Y
  dstat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)
  expect_equal(dstat$discr, 1)

  # small class difference with distance matrix and data produces same result
  sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=1); X <- sim$X; Y <- sim$Y
  dstat.dat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)

  D <- discr.distance(X)
  dstat.dist <- discr.stat(D, Y, is.dist = TRUE)
  expect_equal(dstat.dat$discr, dstat.dist$discr)
})

test_that("One Sample Test is Valid", {
  n = 100; d=5; nsim=50; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(mclapply(1:nsim, function(i) {
    # no true class difference
    set.seed(seed.idx[i])
    sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=0); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(discr.test.one_sample(X, Y, nperm=100)$p.value < alpha)
  }, mc.cores=parallel::detectCores() - 1), use.names=FALSE)
  # check power is near alpha
  expect_lt(abs(mean(res) - alpha), 0.1)
})

test_that("One Sample Test Detects Relationship", {
  n = 100; d=5; nsim=5; alpha=0.1
  set.seed(12345)
  seed.idx <- floor(runif(nsim, 1, 10000))
  res <- unlist(mclapply(1:nsim, function(i) {
    # substantial true class difference
    set.seed(seed.idx[i])
    sim <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=3); X <- sim$X; Y <- sim$Y
    set.seed(seed.idx[i])
    return(discr.test.one_sample(X, Y, nperm=100)$p.value < alpha)
  }, mc.cores=parallel::detectCores() - 1), use.names=FALSE)
  # check power is near 1
  expect_lt(abs(mean(res) - 1), 0.1)
})

test_that("Two Sample Test is Valid", {
  n = 100; d=5; nsim=50; alpha=0.1
  res <- unlist(mclapply(1:nsim, function(i) {
    # no true class difference
    s.g1 <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=0)
    s.g2 <- discr.sims.linear(n=n*3, d=d, K=2, signal.lshift=0)
    g2.out <- list(X=NULL, Y=NULL)
    for (y in unique(s.g1$Y)) {
      idx.g2 <- which(s.g2$Y == y)
      n.y <- sum(s.g1$Y == y)
      g2.out$X <- rbind(g2.out$X, s.g2$X[idx.g2[1:n.y],])
      g2.out$Y <- c(g2.out$Y, s.g2$Y[idx.g2[1:n.y]])
    }
    ord.g1 <- order(s.g1$Y)
    return(discr.test.two_sample(s.g1$X[ord.g1,], g2.out$X, s.g1$Y[ord.g1])$p.value < alpha)
  }, mc.cores=parallel::detectCores() - 1), use.names=FALSE)
  # check power is near alpha
  expect_lt(abs(mean(res) - alpha), 0.1)
})

test_that("Two Sample Test Detects Relationship", {
  n = 100; d=3; nsim=5; alpha=0.1
  res <- unlist(mclapply(1:nsim, function(i) {
    # no true class difference
    s.g1 <- discr.sims.linear(n=n, d=d, K=2, signal.lshift=2, signal.scale=1, non.scale=1)
    s.g2 <- discr.sims.linear(n=n*3, d=d, K=2, signal.lshift=2, signal.scale=2, non.scale=2)
    g2.out <- list(X=NULL, Y=NULL)
    for (y in unique(s.g1$Y)) {
      idx.g2 <- which(s.g2$Y == y)
      n.y <- sum(s.g1$Y == y)
      g2.out$X <- rbind(g2.out$X, s.g2$X[idx.g2[1:n.y],])
      g2.out$Y <- c(g2.out$Y, s.g2$Y[idx.g2[1:n.y]])
    }
    ord.g1 <- order(s.g1$Y)
    return(discr.test.two_sample(s.g1$X[ord.g1,], g2.out$X, s.g1$Y[ord.g1])$p.value < alpha)
  }, mc.cores=parallel::detectCores() - 1), use.names=FALSE)
  # check power is near alpha
  expect_lt(abs(mean(res) - 1), 0.1)
})
