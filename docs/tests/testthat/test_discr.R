library(mgc)
context('Tests for discriminability-related code')

test_that("Test that data ranking works appropriately") {
  n <- 5
  d <- 4
  dat <- array(rnorm(5*4), dim=c(n, d))
  rank_dat <- discr.rank_data(dat)

  # loop over the ranked data, checking that the actual
  # data is ordered in the appropriate way
  for (i in 1:dim(dat)[1]) {
    expect_equal(rank(dat[i,]), rank_dat[i,])
  }
}

test_that("Test that pairwise distances are computed appropriately") {
  n <- 3
  d <- 2
  dat <- array(NaN, dim=c(n, d))
  dat[1,] <- c(1,1)
  dat[2,] <- c(0,0)
  dat[3,] <- c(2,2)

  D <- discr.distance(dat)

  # entries of the distance matrix
  expect_equal(D[1,2], sqrt(2))
  expect_equal(D[2,3], sqrt(2))
  expect_equal(D[1,3], sqrt(4))
  # diagonals are all 0 distance
  for (i in 1:n) {
    expect_equal(D[i,i], 0)
  }
  # symmetry of the distance matrix
  expect_true(all.equal(D, t(D)))
}

# #-----Basic Tests----------------------------------------------#
d <- 16
n <- 2  # number of subjects
nscan <- 2  # number of sessions
X <- array(NaN, dim=c(n*nscan, d))
# two random graphss
x1 <- runif(d)
x2 <- runif(d)
X[1,] <- x1
X[2,] <- x2
X[3,] <- x1
X[4,] <- x2
test_that("Discriminability in perfect match case gets score of 1", {
  labels <- c(0, 1, 0, 1)
  expect_equal(discr.stat(discr.distance(X), ids=labels), 1)
})

test_that("Discriminability in no match case gets score of 1/(nscan*n)", {
  labels <- c(1, 1, 0, 0)
  expect_equal(discr.stat(discr.distance(X), ids=labels), .25)
})
