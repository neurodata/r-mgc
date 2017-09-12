#' KDE of Subject Relationships
#'
#' A function that  returns density estimates of intra and inter subject distances.
#'
#' @param D [n, n]: a distance matrix corresponding to the distances between n scans.
#' @param id [n]: a vector where each entry is the subject id of each of the n scans.
#' @return f_intra [x]: a density estimate of the intra subject relationship.
#' @return f_inter [x]: a density estimate of the inter subject relationship.
#' @return H [1]: the hellinger distance between the intra and inter subject relationships.
#' @author Eric Bridgeford
#' @export
discr.kde_subject <- function(D, id) {
  print("computing hellinger distances and densities...")
  intra <- vector() # array for intra
  inter <- vector() # array for inter
  numsubs <- length(id)
  for (i in 1:(numsubs-1)) {
    for (j in (i+1):numsubs) {
      if (isTRUE(all.equal(id[i],id[j]))) {
        intra <- c(intra, D[i,j]) # if we are intra sub, then append here
      } else {
        inter <- c(inter, D[i,j]) # else we are inter subject
      }
    }
  }
  intert <- density(as.numeric(inter)) # temporary computation so that we can see what range is optimal
  intrat <- density(as.numeric(intra))
  bw_selection <- mean(c(intert$bw, intrat$bw))
  minx <- min(c(intrat$x, intert$x)) #compute min
  maxx <- max(c(intrat$x, intert$x)) #and max vals
  f_intra <- density(intra, from=minx, to=maxx, bw = bw_selection) # compute the density for the predetermined
  f_inter <- density(inter, from=minx, to=maxx, bw = bw_selection) # range(from, to)

  H = norm(sqrt(f_intra$y)-sqrt(f_inter$y), type="2")/sqrt(2) # hellinger distance from the yvals

  return(list(f_intra, f_inter, H))
}
