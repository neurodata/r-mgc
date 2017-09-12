#' Distance
#'
#' A function that returns a distance matrix given a collection of graphs. Note that this function uses the frobenius norm to compute the distance between all pairs of graphs.
#'
#' @import stats
#' @param graphs: [[n]][nroi, nroi] a collection of graphs for n subjects, where each graph is a square matrix.
#' @param method='F': the method for computing distances.
#' \describe{
#'    \item{'F'}{use the frobenius norm between pairs of graphs}
#'    \item{'H'}{use the hellinger distance between pairs of graphs}
#' }
#' @return dist [n, n]: a matrix indicating the pairwise distances between all graphs passed in.
#' @author Eric Bridgeford
#' @export
discr.distance <- function(graphs, method='F') {
  if (method == 'F') {
    graphs <- fmriu.list2array(graphs)
    dim_graphs <- dim(graphs) # get the dims of the graphs
    # expect dim_graphs[1] to be nrois, dim_graphs[2] to be nois, dim_graphs[3] to be numsubs
    reshape_graphs <- t(array(graphs, dim=c(dim_graphs[1]*dim_graphs[2], dim_graphs[3])))
    dist_graphs <- dist(reshape_graphs, diag=TRUE, upper=TRUE) # use stats dist function
    D <- array(matrix(as.matrix(dist_graphs)), dim=c(dim_graphs[3],dim_graphs[3]))
  } else if (method == 'H') {
    nsub <- length(div)
    D <- array(NaN, dim=c(nsub, nsub))
    for (sub1 in 1:nsub) {
      for (sub2 in 1:nsub) {
        D[sub1, sub2] <- discr.hell_dist(div[[sub1]], div[[sub2]])
      }
    }
  } else {
    stop('You have not passed a valid method. Please select either [\'F\', \'H\'].')
  }
  return(D)
}
