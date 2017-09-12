#' Graph Discriminability
#'
#' A function that performs discriminability on a set of graphs.
#'
#' @import fmriutils
#' @import Rmisc
#' @import ggplot2
#' @import reshape2
#' @param graphs [[n]][nroi, nroi]: the graphs for each of the n subjects with nroi rois as a list of adjacency matrices.
#' @param ids [n]: the ids for each scan corresponding to the signal from above.
#' @param rank=FALSE : a boolean indicating whether to do unranked (FALSE) or ranked (TRUE).
#' @param font.size=15: the default font size.
#' @param method='F' the method to use for the distance computation.
#' \describe{
#'   \item{'F'}{Frobenius norm between pairs of graphs}
#'   \item{'H'}{Hellinger distance between pairs of graphs}
#' }
#' @return d [1]: the discriminability statistic for the data.
#' @return dist [n, n]: the distance matrix associated with the data.
#' @return distplot : a plot of the distance matrix.
#' @return kdeplot : a plot of the density estimate of the intra vs inter subject distances.
#' @return combinedplot : a multiplot showing the dist plot and the kde plot.
#' @author Eric Bridgeford
#' @export
discr.graphs2discr <- function(graphs, ids, rank=FALSE, font.size=15, method='F') {
  if (rank==TRUE) {
    graphs <- discr.rank_edges(graphs)
  }

  D <- discr.distance(graphs, method=method)
  discrstat <- discr.discr(discr.rdf(D, ids), thresh=.01)

  kdeobj <- discr.kde_subject(D, ids)
  kde_dist <- data.frame(x=kdeobj[[1]]$y, y=kdeobj[[2]]$y, distance=kdeobj[[1]]$x)
  colnames(kde_dist) <- c("intra", "inter", "distance")
  meltkde <- melt(kde_dist, id="distance")
  colnames(meltkde) <- c("distance", "Relationship", "Probability")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  distance_plot <- ggplot(melt(D), aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("darkblue","blue","purple","green","yellow"),
                         name="dist") +
    xlab("Scan") +
    ylab("Scan") +
    ggtitle(sprintf('Distance Matrix, d=%.4f', discrstat)) +
    theme(text=element_text(size=font.size))

  kde_plot <- ggplot() +
    geom_ribbon(data=meltkde, aes(x=distance, ymax=Probability, fill=Relationship), ymin=0, alpha=0.5) +
    ggtitle('Density Estimate') +
    theme(text=element_text(size=font.size))
  dual_plot <- multiplot(distance_plot, kde_plot, layout=matrix(c(1,2), nrow=1, byrow=TRUE))
  return(list(d = discrstat, dist = D, distplot = distance_plot, kdeplot = kde_plot, combinedplot = dual_plot))
}
