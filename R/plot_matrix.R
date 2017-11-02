#' Matrix plot
#'
#' A function that plots a matrix.
#'
#' @import ggplot2
#' @import reshape2
#' @param mtx [nroi, nroi] or [nt, nroi]: the input. Can either be a square matrix of featuresxfeatures, or a data matrix of observationsxfeatures.
#' @param title="": the title for the square plot.
#' @param xlabel="": the x label for the square plot.
#' @param ylabel="": the y label for the square plot.
#' @param legend.name="": the legend title for the square plot.
#' @param legend.show=TRUE: whether to show the legend on the plot.
#' @param font.size=12: the default font size for the plot text. Axis/legend text is font.size - 2.
#' @author Eric Bridgeford
#' @export
mgc.plot.plot_matrix <- function(mtx, title="",xlabel="", ylabel="", legend.name="metric", legend.show=TRUE,
                                 font.size=12, limits=NaN) {
  dm <- reshape2::melt(mtx)
  if (is.nan(limits)) {
    limits <- c(min(mtx), max(mtx))
  }
  colnames(dm) <- c("x", "y", "value")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  sqplot <- ggplot2::ggplot(dm, aes(x=x, y=y, fill=value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours=jet.colors(7), name=legend.name, limits=limits) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(title)
  if (legend.show) {
    sqplot <- sqplot +
      ggplot2::theme(text=element_text(size=font.size))
  } else {
    sqplot <- sqplot +
      ggplot2::theme(text=element_text(size=font.size, legend.position="none"))
  }
  return(sqplot)
}
