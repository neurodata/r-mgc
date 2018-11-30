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
  sim <- discr.sims.linear(n=n, d=d, K=2, mean.scale=10); X <- sim$X; Y <- sim$Y
  dstat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)
  expect_equal(dstat$discr, 1)

  # small class difference with distance matrix and data produces same result
  sim <- discr.sims.linear(n=n, d=d, K=2, mean.scale=1); X <- sim$X; Y <- sim$Y
  dstat.dat <- discr.stat(X, Y, remove.isolates = FALSE, is.dist = FALSE)

  D <- discr.distance(X)
  dstat.dist <- discr.stat(D, Y, is.dist = TRUE)
  expect_equal(dstat.dat$discr, dstat.dist$discr)
})
