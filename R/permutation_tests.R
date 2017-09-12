#' Discriminability One Sample Permutation Test
#'
#' A function that permutes the labels of a distance matrix to obtain an empirical pvalue associated with whether the raw score is due to random chance.
#'
#' @param D [nxn] the distance matrix to run a permutation test for. An [nxn] matrix.
#' @param labels [n] the labels organized appropriately with the distance matrix. Label 1 should correspond to the first column, 2 the second, and so on.
#' @param nrep=100 the number of permutations to perform.
#' @param verbose=FALSE whether to print the itereation numbers.
#' @return srel the relative, unpermuted discriminability you want to see is significant.
#' @return null the discriminability scores of the permuted data.
#' @return pval the pvalue associated with the permutation test.
#' @author Shangsi Wang
#' @export
discr.test.one_sample <- function(D, labels, nperm=100, verbose=FALSE) {
  # test whether discriminability differs from 0.5
  N <- dim(D)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  tr <- discr.(discr.rdf(D, ids))

  nr <- rep(0,nperm)
  for (i in 1:nperm){
    if (verbose) {
      print(i)
    }
    nr[i] <- mean(discr.rdf(D, ids[sample(N)]))
  }
  result <- list()
  result$srel <- tr
  result$null <- sort(nr)
  result$pval <- sum(nr>tr)/nperm
  return(result)
}

#' Discriminability Two Sample Permutation Test
#'
#' A function that permutes the labels of a distance matrix to obtain an empirical pvalue associated with whether the raw score is due to random chance.
#'
#' @param D1 [nxn] the first distance matrix to run a permutation test for. An [nxn] matrix.
#' @param D2 [nxn] the second distance matrix to run a permutation test for. An [nxn] matrix.
#' @param labels [n] the labels organized appropriately with the distance matrices. Label 1 should correspond to the first column, 2 the second, and so on.
#' @param verbose=FALSE whether to print the itereation numbers.
#' @param nperm=100 the number of permutations to perform.
#' @return pval the pvalue associated with the permutation test.
#' @author Shangsi Wang
#' @export
discr.test.two_sample <- function(dist1, dist2, ids, nperm=100, verbose=FALSE){
  # test two discriminability are the same
  N1 <- dim(dist1)[1]
  N2 <- dim(dist2)[1]

  if (is.null(N1) || is.null(N2) ) {
    stop('Invalid datatype for dist1 or dist2')
  }
  if (N1 != N2) {
    stop('The dimension of dist1 and dist2 do not match')
  }

  disct1 <- matrix(0,N1,2)
  disct2 <- matrix(0,N1,2)

  for (i in 1:N1){
    disct1[i,] <- dis_vec(dist1[i,],i,ids)
    disct2[i,] <- dis_vec(dist2[i,],i,ids)
  }

  disct1 <- disct1[!is.na(disct1[,1]),]
  disct2 <- disct2[!is.na(disct2[,1]),]
  N1 <- dim(disct1)[1]
  N2 <- dim(disct2)[1]

  if (N1 != N2) {
    stop('The dimension of dist1 and dist2 do not match')
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
  return (list(pval=pvalue))

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
    stop('Invalid datatype for dist1 or dist2')
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



