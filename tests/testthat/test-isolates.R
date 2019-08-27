context("isolates")

test_that("Removing Isolates - d == 1, 2 Class", {
  X <- as.matrix(c(0, 1, 2, 3, 4)); Y <- c(1, 1, 2, 1, 1)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, as.matrix(c(0, 1, 3, 4)))
  expect_equal(purged$Y, c(1, 1, 1, 1))

  Y <- c(2, 1, 1, 1, 1)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, as.matrix(c(1, 2, 3, 4)))
  expect_equal(purged$Y, c(1, 1, 1, 1))

  Y <- c(1, 1, 1, 1, 2)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, as.matrix(c(0, 1, 2, 3)))
  expect_equal(purged$Y, c(1, 1, 1, 1))
})

test_that("Removing Isolates - d > 1, 2 Class", {
  n = 5; d = 10
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X; Y <- c(2, 1, 1, 1, 1)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,])
  expect_equal(purged$Y, c(1, 1, 1, 1))
  Y <- c(1, 1, 2, 1, 1)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,])
  expect_equal(purged$Y, c(1, 1, 1, 1))

  Y <- c(1, 1, 1, 1, 2)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,])
  expect_equal(purged$Y, c(1, 1, 1, 1))
})

test_that("Removing 2 Isolated Classes - d == 1, 3 Class", {
  n = 6; d = 1
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X; Y <- c(2, 1, 1, 1, 1, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])

  Y <- c(1, 1, 2, 1, 1, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])

  Y <- c(1, 1, 1, 1, 2, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])
})

test_that("Removing 2 Isolated Classes - d > 1, 3 Class", {
  n = 6; d = 10
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X; Y <- c(2, 1, 1, 1, 1, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated sample is removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])

  Y <- c(1, 1, 2, 1, 1, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated samples removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])

  Y <- c(1, 1, 1, 1, 2, 3)
  purged <- remove.isolates(X, Y)
  # validate that isolated samples removed correctly
  expect_equal(purged$X, X[Y == 1,,drop=FALSE])
  expect_equal(purged$Y, Y[Y == 1])
})

test_that("All data is isolated", {
  n = 6; d = 10
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X; Y <- seq(1, n)
  purged <- remove.isolates(X, Y)
  # validate that isolated samples removed correctly
  expect_equal(purged$X, matrix(0, nrow=0, ncol=d))
  expect_equal(purged$Y, integer(0))
})

test_that("Distance Matrix Isolate Removal", {
  n = 6; d = 10
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X;
  D <- mgc.distance(X); Y <- c(2, 1, 1, 1, 1, 1)
  purged <- remove.isolates(D, Y, is.dist=TRUE)
  # validate that isolated samples removed correctly
  expect_equal(as.numeric(purged$X), as.numeric(D[which(Y == 1), which(Y == 1)]))
  expect_equal(purged$Y, c(1, 1, 1, 1, 1))

  n = 7; d = 10
  sim <- discr.sims.linear(n=n, d=d, K=2); X <- sim$X;
  D <- mgc.distance(X); Y <- c(2, 1, 1, 1, 3, 1, 1)
  purged <- remove.isolates(D, Y, is.dist=TRUE)
  # validate that isolated samples removed correctly
  expect_equal(as.numeric(purged$X), as.numeric(D[which(Y == 1), which(Y == 1)]))
  expect_equal(purged$Y, c(1, 1, 1, 1, 1))

})
