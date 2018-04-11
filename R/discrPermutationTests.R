#' Discriminability One Sample Permutation Test
#'
#' A function that permutes the labels of a distance matrix to obtain an empirical pvalue associated with whether the raw score is due to random chance.
#'
#' @param D \code{[n x n]} the distance matrix to run a permutation test for, for \code{n} samples.
#' @param ids \code{n} the labels of each of the \code{n} samples, in the same ordering as elements of the distance matrix. Label 1 should correspond to the first column, 2 the second, and so on.
#' @param nperm the number of permutations to perform. Defaults to \code{100}.
#' @param verbose whether to print the itereation numbers. Defaults to \code{FALSE}.
#' @return A list containing the following:
#' \item{\code{srel}}{the relative, unpermuted discriminability you want to see is significant.}
#' \item{\code{null}}{the discriminability scores of the permuted data.}
#' \item{\code{pval}}{the pvalue associated with the permutation test.}
#' @author Shangsi Wang and Eric Bridgeford
#' @export
discr.test.one_sample <- function(D, ids, nperm=100, verbose=FALSE) {
  # test whether discriminability differs from 0.5
  N <- dim(D)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  tr <- discr.mnr(discr.rdf(D, ids))

  nr <- rep(0,nperm)
  for (i in 1:nperm){
    if (verbose) {
      print(i)
    }
    nr[i] <- discr.mnr(discr.rdf(D, ids[sample(N)]))
  }
  result <- list()
  result$srel <- tr
  result$null <- sort(nr)
  result$pval <- sum(nr>tr)/nperm
  return(result)
}

#' Discriminability Two Sample Permutation Test
#'
#' A function that takes two distance matrices and produces a p-value associated with whether or not the distance matrices differ significantly.
#'
#' @param D1 \code{[n x n]} the first distance matrix to run a permutation test for, for \code{n} samples.
#' @param D2 \code{[n x n]} the second distance matrix to run a permutation test for, for \code{n} samples.
#' @param ids \code{n} the labels of each of the \code{n} samples, in the same ordering as elements of the distance matrix. Label 1 should correspond to the first column, 2 the second, and so on.
#' @param nperm the number of permutations to perform. Defaults to \code{100}.
#' @param verbose whether to print the itereation numbers. Defaults to \code{FALSE}.
#' @return the pvalue associated with the permutation test.
#' @author Shangsi Wang and Eric Bridgeford
#' @export
discr.test.two_sample <- function(D1, D2, ids, nperm=100, verbose=FALSE){
  # test two discriminability are the same
  N1 <- dim(D1)[1]
  N2 <- dim(D2)[1]

  if (is.null(N1) || is.null(N2) ) {
    stop('Invalid datatype for D1 or D2')
  }
  if (N1 != N2) {
    stop('The dimension of D1 and D2 do not match')
  }

  disct1 <- matrix(0,N1,2)
  disct2 <- matrix(0,N1,2)

  for (i in 1:N1){
    disct1[i,] <- discr.test.dis_vec(D1[i,],i,ids)
    disct2[i,] <- discr.test.dis_vec(D2[i,],i,ids)
  }

  disct1 <- disct1[!is.na(disct1[,1]),]
  disct2 <- disct2[!is.na(disct2[,1]),]
  N1 <- dim(disct1)[1]
  N2 <- dim(disct2)[1]

  if (N1 != N2) {
    stop('The dimensions of D1 and D2 do not match')
  }

  tcount <- sum(disct1[,2])
  tdif <- (sum(disct1[,1] * disct1[,2]) - sum(disct2[,1] * disct2[,2])) / tcount

  ndif <- rep(0,nperm)
  ndisct1 <- matrix(0,N1,2)
  ndisct2 <- matrix(0,N1,2)
  for (i in 1:nperm){
    if (verbose) {
      print(i)
    }
    ind <- rbinom(N1,1,0.5) == 1
    for (j in 1:N1){
      ndisct1[ind,] <- disct1[ind,]
      ndisct1[!ind,] <- disct2[!ind,]
      ndisct2[ind,] <- disct2[ind,]
      ndisct2[!ind,] <- disct1[!ind,]
    }
    ndif[i] <- (sum(ndisct1[,1] * ndisct1[,2]) - sum(ndisct2[,1] * ndisct2[,2])) / tcount
  }
  pvalue <- (sum(ndif > abs(tdif)) + 0.5 * sum(ndif == abs(tdif))) / nperm
  return (pval=pvalue)
}

discr.test.dis_vec<-function(distvec,i,ids){
  N <- length(distvec)
  ind <- which(grepl(ids[i],ids))
  rdf <- c()
  count <- 0
  for (j in ind) {
    if (j != i) {
      di <- distvec
      d <- di[j]
      di[ind] <- Inf
      count <- count + 1
      rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind))
    }
  }
  return(c(mean(rdf),count))
}

discr.test.reliability_bootstrap<-function(ids1,obs1,nrep=100){
  # test two discriminability are the same
  N <- length(ids1)
  if (is.null(N) ) {
    stop('Invalid datatype for D1 or D2')
  }

  rels <- rep(0,nrep)
  uids1 <- unique(ids1)
  bobs <- obs1
  pduct <- obs1 %*% t(obs1)



  for (i in 1:nrep){
    bids1 <- c()
    cind <- 1

    ind1v<-c()
    ind2v<-c()
    tv<-c()

    for (m in 1:(length(uids1))){
      t <- runif(1)
      nmeas <- length(which(grepl(uids1[m],ids1)))
      nmeas1<-0
      nmeas2<-0
      while(min(nmeas1,nmeas2)<nmeas){
        twoids <- sample(uids1,2)
        nmeas1 <- length(which(grepl(twoids[1],ids1)))
        nmeas2 <- length(which(grepl(twoids[2],ids1)))
      }
      ind1<-sample(which(grepl(twoids[1],ids1)),nmeas)
      ind2<-sample(which(grepl(twoids[2],ids1)),nmeas)
      ind1v<-c(ind1v,ind1)
      ind2v<-c(ind2v,ind2)
      tv <- c(tv,rep(t,nmeas))

      nobs <- obs1[ind1,]*t + obs1[ind2,] * (1-t)
      bobs[cind:(cind+nmeas-1),] <- nobs
      cind <- cind + nmeas
      bids1 <- c(bids1, rep(paste('pseudo',m+10000,sep=""),nmeas))
    }


    distmat1 <- matrix(0,N,N)



    if (ncol(obs1) > 1){
      for (j in 1:N){
        for(k in j:N){
          dv <- rep(0,N)
          dv[ind1v[j]] <- dv[ind1v[j]] + tv[j]
          dv[ind2v[j]] <- dv[ind2v[j]] + (1-tv[j])
          dv[ind1v[k]] <- dv[ind1v[k]] - tv[k]
          dv[ind2v[k]] <- dv[ind2v[k]] - (1-tv[k])
          distmat1[j,k] <- t(dv) %*% pduct %*% dv
        }
      }
      distmat1 <- sqrt(distmat1)
      distmat1<-distmat1 + t(distmat1)
    } else {
      bO <- obs1[ind1v,] * tv +  obs1[ind2v,] *(1-tv)
      distmat1 <- as.matrix(dist(bO,diag = T, upper = T))
      distmat1 <- distmat1+ t(distmat1)
    }


    rels[i] <- mean(discr.rdf(distmat1,bids1))
  }
  return(rels)
}
