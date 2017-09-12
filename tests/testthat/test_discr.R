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

# #-----Basic Tests----------------------------------------------#
# nroi <- 2
# n <- 2
# nscan <- 2
# graphs <- array(NaN, dim=c(nroi, nroi, n*nscan))
# # two random graphss
# g1 <- array(runif(nroi*nroi), dim=c(nroi, nroi))
# g2 <- array(runif(nroi*nroi), dim=c(nroi, nroi))
# graphs[,,1] <- g1
# graphs[,,2] <- g2
# graphs[,,3] <- g1
# graphs[,,4] <- g2
# graphs <- fmriu.array2list(graphs)
# test_that("Discriminability in perfect match case gets score of 1", {
#   labels <- c(0, 1, 0, 1)
#   expect_equal(discr.discr(discr.rdf(discr.distance(graphs), labels)), 1)
# })
#
# test_that("Discriminability in no match case gets score of 1/(nscan*n)", {
#   labels <- c(1, 1, 0, 0)
#   expect_equal(discr.discr(discr.rdf(discr.distance(graphs), labels)), .25)
# })
#
# test_that("Discriminability via graphs2discr driver works in perfect match case", {
#   labels <- c(0, 1, 0, 1)
#   expect_equal(discr.graphs2discr(graphs, labels, rank=FALSE), 1)
# })
#
#
#
# test_that("Discriminability via graphs2discr driver works in no match case gets score of 1/(nscan*n)", {
#   labels <- c(0, 0, 1, 1)
#   expect_equal(discr.graphs2discr(graphs, labels, rank=FALSE), .25)
# })
