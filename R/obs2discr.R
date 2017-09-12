#' Frequency Discriminability
#'
#' A function that performs discriminability in the frequency domain.
#'
#' @import fmriutils
#' @import ggplot2
#' @import reshape2
#' @import Rmisc
#' @param signal [[n]][nt, nroi] the signal for each of the n subjects, containing an array of nt observations for nroi rois.
#' @param ids [n] the ids for each scan corresponding to the signal from above.
#' @param tr=NaN [1] the repetition time of the dataset. NaN for none.
#' @param lc=NaN [1] the lower cutoff for highpass filtering. NaN for none.
#' @param hc=NaN [1] the high cutoff for low-pass filtering. NaN for none.
#' @param spec='amp' the spectrum to work with. 'amp' for amplitude, 'pow' for power.
#' @param rank=FALSE a boolean indicating whether to do unranked (FALSE) or ranked (TRUE) graphs, whereby all the edge-weights gare ranked before computing distances.
#' @param font.size=10 the default font size for the plot text.
#' @param method='H' the method to use for the distance computation.
#' \describe{
#'   \item{'F'}{Frobenius norm between pairs of graphs}
#'   \item{'H'}{Hellinger distance between pairs of graphs}
#' }
#' @return d [1] the discriminability statistic for the data.
#' @return dist [n, n] the distance matrix associated with the data.
#' @return distplot a plot of the distance matrix.
#' @return kdeplot a plot of the density estimate of the intra vs inter subject distances.
#' @return combinedplot a multiplot showing the dist plot and the kde plot.
#' @author Eric Bridgeford
#' @export
discr.freq.obs2discr <- function(signal, ids, tr, lc=NaN, hc=NaN, spectrum='amp', rank=FALSE, font.size=15, method='H') {
  spec_sig <- fmriu.freq.obs2div(signal, tr=tr, lc=lc, hc=hc, spectrum=spec)
  div <- freq2div(spec_sig)
  if (rank==TRUE) {
    div <- discr.rank_edges(div)
  }

  D <- discr.distance(div, method)
  discrspec <- discr.discr(discr.rdf(D, ids), thresh=0.01)

  kdeobj <- discr.kde_subject(D, ids)
  kde_dist <- data.frame(x=kdeobj[[1]]$y, y=kdeobj[[2]]$y, distance=kdeobj[[1]]$x)
  colnames(kde_dist) <- c("intra", "inter", "distance")
  meltkde <- melt(kde_dist, id="distance")
  colnames(meltkde) <- c("distance", "Relationship", "Probability")

  distance_plot <- ggplot(melt(D), aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("darkblue","blue","purple","green","yellow"),name="distance") +
    xlab("Scan") +
    ylab("Scan") +
    ggtitle(sprintf('Distance Matrix between Pairs of Scans, d=%.4f', discrspec)) +
    theme(text=element_text(size=font.size))
  kde_plot <- ggplot() +
    geom_ribbon(data=meltkde, aes(x=distance, ymax=Probability, fill=Relationship), ymin=0, alpha=0.5) +
    ggtitle('Density Estimate of Scan Distances') +
    theme(text=element_text(size=font.size))
  dual_plot <- multiplot(distance_plot, kde_plot, layout=matrix(c(1,2), nrow=1, byrow=TRUE))
  return(list(d=discrspec, dist=D, distplot=distance_plot, kdeplot = kde_plot, combinedplot=dual_plot))
}

#' Time Domain Discriminability
#'
#' A function that performs discriminability in the time domain using the correlation.
#'
#' @import fmriutils
#' @import Rmisc
#' @import ggplot2
#' @import reshape2
#' @param signal [[n]][nt, nroi]: the signal for each of the n subjects, containing an array of nt observations,  for nroi rois. Alternatively, can be the graphs for each of the subjects.
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
discr.time.obs2discr <- function(signal, ids, rank=FALSE, graphs=FALSE, font.size=15, method='F') {
  corr <- fmriu.time.obs2corr(signal)
  if (rank==TRUE) {
    corr <- discr.rank_edges(corr)
  }

  D <- discr.distance(corr, method=method)
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
