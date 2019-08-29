context("distance")
# testing for distance function

test_that("Distance Function Simple Example - Orthogonal 2 samples", {
  X <- rbind(c(0, 1), c(1, 0))
  D <- mgc.distance(X)
  expect_equal(as.numeric(D), as.numeric(cbind(c(0, sqrt(2)), c(sqrt(2), 0))))
})

test_that("Distance Function Simple Example - Same 2 samples", {
  X <- rbind(c(1, 0), c(1, 0))
  D <- mgc.distance(X)
  expect_equal(as.numeric(D), as.numeric(cbind(c(0, 0), c(0, 0))))
})
