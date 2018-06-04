#' Multiscale Network Test
#'
#' This function provides multiscale test statistics and p-values, and those under optimal signal
#' in testing independence between edge connectivity in graph \code{G} and nodal attributes \code{X} using diffusion map embeddings. Default test method is \code{mgc}.
#'
#'
#' @param G the input igraph object.
#' @param X the input nodal attributes.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{uses the MGC global correlation.}
#'    \item{\code{'dcor'}}{uses the dcor global correlation.}
#'    \item{\code{'hhg'}}{uses Heller-Heller-Gorfine (HHG) tests.}
#' }
#' @param t.max the maximum diffusion time. Diffusion time starts from 0,1,.. \code{t.max}. Defaults to \code{t.max} = 5.
#' @param n.perm the number of permutation samples.  Defaults to \code{n.perm} = 100.
#' @param default.q dimension of graph embeddings. Defaults to \code{default.q} = 10.
#' @param default.t default diffusion time at optimal when no significant signal was found. Defaults to \code{default.t} = 3.
#'
#' @return A list containing the following:
#' \item{\code{pval}}{is p-value under diffusion map at \code{optimal.t}.}
#' \item{\code{optimal.t}}{diffusion time at smoothed maximum among multiscale statistics.}
#' \item{\code{dimension}}{a dimension of diffusion map embedding.}
#' \item{\code{localmap}}{a list of local correlation map (in case of \code{option} = '\code{mgc}').}
#' \item{\code{localscale}}{a list of optimal local scale (row, column) in each local correlation map (in case of \code{option} = '\code{mgc}').}
#'
#' @author Youjin Lee
#'
#' @export
#'
#' @examples
#' library(mgc)
#' library(igraphdata)
#' data(karate)
#' G <- karate
#' X <- V(karate)$Faction
#' result <- NetworkTest(G, X, option = 'mgc', t.max = 5)
#' print(result$pval) # p-value
#' print(result$optimal.t) # optimal diffusion time
#'
NetworkTest = function(G, X, option = 'mgc', t.max = 5, n.perm = 100, default.q = 10, default.t = 3){

  test.result = DiffusionTest(G, X, option, t.max, n.perm, default.q)

  if(length(test.result$stats) < 3){
    optimal.t = 1;
    optimal.pval = mean(test.result$stats.null[,optimal.t] > test.result$stats[optimal.t])
  }
  optimal.result = print.stat.optimal(list(test.result$stats, test.result$stats.null), default.t+1)
  optimal.t = optimal.result$alt.t - 1
  optimal.pval = optimal.result$pvalue

  return(list(pval = optimal.pval, optimal.t = optimal.t, dimension = test.result$dimension,
              localmap = test.result$localmap, localscale = test.result$localscale))

}

###################
#' Distance-based graph test using diffusion map embedding
#'
#'
#' An auxiliary function to print out multiscale statistics and those under the permutation sample
#' along with local correlations map across diffusion time.
#'
#' @importFrom  meda getElbows
#' @importFrom HHG hhg.test
#' @importFrom energy dcor.ttest
#'
#' @param G the input igraph object.
#' @param X the input nodal attributes.
#' @param option is a string that specifies which global correlation to build up-on. Defaults to \code{'mgc'}.
#' \describe{
#'    \item{\code{'mgc'}}{uses the MGC global correlation.}
#'    \item{\code{'dcor'}}{uses the dcor global correlation.}
#'    \item{\code{'hhg'}}{uses Heller-Heller-Gorfine (HHG) tests.}
#' }
#' @param t.max the maximum diffusion time. It starts from 0, and default is \code{t.max} = 5.
#' @param n.perm the number of permutation samples.  Defaults to \code{n.perm} = 100.
#' @param default.q a dimension of graph embeddings. Defaults to \code{default.q} = 10.
#'
#' @return A list containing the following:
#' \item{\code{stats}}{a set of multiscale statistics across diffusion times.}
#' \item{\code{stats.null}}{a \code{[n.perm x (t.max + 1)]} matrix of multiscale statistics under the null.}
#' \item{\code{dimension}}{a dimension of diffusion map embedding.}
#' \item{\code{localmap}}{a list of local correlation map (in case of \code{option} = 'mgc').}
#' \item{\code{localscale}}{a list of optimal local scale (row, column) in each local correlation map (in case of \code{option} = 'mgc').}
#'
#'
#' @author Youjin Lee
#'
#' @export
#'
DiffusionTest = function(G, X, option = 'mgc', t.max = 5, n.perm = 100, default.q = 10){

  A = as.matrix(get.adjacency(G))
  D = diag(pmin( (rowSums(A))^(-1/2) , 1))
  P = D %*% A %*% D
  t.range = c(0:t.max)

  mgc.stat = c(); dcor.stat = c(); hhg.stat = c()
  localmap = list(); localscale = list()

  diffusion.q = try(min(max(getElbows(abs(print.lambda(P, times = 1)[[1]]), plot = FALSE, n = 2)), nrow(A)-1),silent = TRUE)
  if(class(diffusion.q) == "try-error") diffusion.q = default.q

  for(s in 1:length(t.range)){
    U  =  dmap.q(P, t.range[s], diffusion.q)[[1]]
    Dx = as.matrix(dist((U)), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)

    if(option == 'mgc'){
      mgc = mgc.sample(Dx, Dy, option = 'mgc')
      mgc.stat[s] = mgc$statMGC
      localmap[[s]] = mgc$localCorr
      localscale[[s]] = c(mgc$optimalScale$x, mgc$optimalScale$y)
    }

    if(option == 'dcor'){
      dcor = dcor.ttest(Dx, Dy, distance = TRUE)
      dcor.stat[s] = dcor$statistic
    }


    if(option == 'hhg'){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      hhg.stat[s] = hhg$sum.lr
    }
  }

  # permutation
  tmp = matrix(0, nrow = n.perm, ncol = length(t.range))
  for(r in 1:n.perm){
    for(s in 1:length(t.range)){
      per=sample(length(X));
      newX = X[per]
      Dx = as.matrix(dist(U), diag = TRUE, upper = TRUE)
      Dy = as.matrix(dist(newX), diag = TRUE, upper = TRUE)
      if(option == 'mgc'){
        mgc = mgc.sample(Dx, Dy, option = 'mgc')
        tmp[r,s] = mgc$statMGC
      }
      if(option == 'dcor'){
        dcor = dcor.ttest(Dx, Dy, distance = TRUE)
        tmp[r,s] = dcor$statistic
      }
      if(option == 'hhg'){
        hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
        tmp[r,s] = hhg$sum.lr
      }
    }
  }

  if(option == 'mgc') return(list(stats = mgc.stat, stats.null = tmp, dimension = diffusion.q,
                                  localmap = localmap, localscale = localscale))
  if(option == 'dcor') return(list(stats = dcor.stat, stats.null = tmp, dimension = diffusion.q,
                                   localmap = localmap, localscale = localscale))
  if(option == 'hhg') return(list(stats = hhg.stat, stats.null = tmp, dimension = diffusion.q,
                                  localmap = localmap, localscale = localscale))
}


######################
#' Eigenvalues of diffusion map
#'
#' This auxiliary function provides an eigenvalues of the diffusion maps at given diffusion time.
#'
#' @usage
#' print.lambda(L, times = 1)
#'
#' @importFrom DTMCPack statdistr
#' @importFrom Matrix Diagonal
#'
#' @param L [n x n] Laplacian matrix.
#' @param times a vector of diffusion times. Defaults to 1.
#'
#' @return a list of eigenvalues at each diffusion time.
#' @export
#'
#' @author Youjin Lee
#'
#'
print.lambda = function(L, times = 1){
  maps = list()
  pi = statdistr(L)
  isolated = which(pi <= 0)
  pi[isolated] = 0
  pi.mat1 = as.matrix(Diagonal(length(pi), pi^(1/2)))
  pi.mat2 = as.matrix(Diagonal(length(pi), pi^(-1/2)))

  pi.mat1[,isolated] = 0
  pi.mat2[,isolated] = 0

  # symmetric kernel
  Q = pi.mat1 %*% L %*% pi.mat2

  Q[,isolated] = 0
  Q[,isolated] = 0

  lambda = Re(eigen(Q)$values)[Im(eigen(Q)$value) == 0]
  psi = Re(eigen(Q)$vectors)[,Im(eigen(Q)$value) == 0]

  phi = pi.mat2 %*% psi

  Lambda = as.matrix(Diagonal(length(lambda), lambda))

  lambda.set = list()
  for(t in 1:length(times)){
    lambda.set[[t]] = lambda^(times[t])
  }
  return(lambda.set)
}

#########

#' Diffusion maps
#'
#' This function provides a set of diffusion maps.
#'
#' @usage
#' dmap.q(L, times = 1, q = 10)
#'
#' @importFrom DTMCPack statdistr
#' @importFrom Matrix Diagonal
#'
#' @param L [n x n] Laplacian matrix.
#' @param times a vector of diffusion times. Defaults to 1.
#' @param q dimension of diffusion map embeddings.
#'
#' @return a list of diffusion map embeddings with length [q].
#' @export
#'
#' @author Youjin Lee
#'
#' @examples
#'
#' library(igraphdata)
#' library(mgc)
#' data(karate)
#' G <- karate
#' A <- as.matrix(get.adjacency(G))
#' D <- diag(pmin( (rowSums(A))^(-1/2) , 1))
#' P <- D %*% A %*% D
#' diffusion.map <- dmap.q(P, c(0:5), q = 5)
#'
dmap.q = function(L, times = 1, q = 10){

  maps = list(); maps.q = list()
  pi = statdistr(L) # stationary distribution (probability)
  isolated = which(pi <= 0)
  pi[isolated] = 0
  pi.mat1 = as.matrix(Diagonal(length(pi), pi^(1/2)))
  pi.mat2 = as.matrix(Diagonal(length(pi), pi^(-1/2)))

  pi.mat1[,isolated] = 0
  pi.mat2[,isolated] = 0
  # symmetric kernel
  Q = pi.mat1 %*% L %*% pi.mat2

  Q[,isolated] = 0
  Q[,isolated] = 0

  lambda = Re(eigen(Q)$values)[Im(eigen(Q)$value) == 0]
  psi = Re(eigen(Q)$vectors)[,Im(eigen(Q)$value) == 0]


  phi = pi.mat2 %*% psi

  Lambda = as.matrix(Diagonal(length(lambda), lambda))

  # each row of maps is a diffusion maps of each vertex
  for(t in 1:length(times)){
    maps[[t]] <- phi %*% Lambda^(times[t])
    maps.q[[t]] <- maps[[t]][,1:q]
  }

  return(maps.q)
}


#' Smoothed maximum statistics
#'
#'
#' This function prints out the statistic at smoothed maximum diffusion time,
#'  its p-value and the statistics under the null based on the permutation.
#'
#' @usage
#' print.stat.optimal(stat.list, default.t = 3)
#'
#' @param stat.list a list containing following statistics:
#' \describe{
#' \item{\code{stats}}{\code{[t]} statistics across diffusion times.}
#' \item{\code{stats.null}}{\code{[n.perm, t]} statistics under the permutation sample.}
#' }
#' @param default.t a default diffusion time when no significant signal was found. Defaults to \code{default.t} = 3.
#'
#' @return A list containing the following:
#' \describe{
#' \item{\code{pval}}{is p-value under diffusion map at \code{alt.t}.}
#' \item{\code{alt.stat}}{a test statistic at \code{alt.t}.}
#' \item{\code{alt.t}}{diffusion time at smoothed maximum among multiscale statistics.}
#' \item{\code{null.stat}}{test statistics at \code{null.t} upon permutation sample.}
#' \item{\code{null.t}}{diffusion time at smoothed maximum among multiscale statistics under the null.}
#' }
#'
#' @export
#'
#' @author Youjin Lee
#'
#'
#'
#'
print.stat.optimal = function(stat.list, default.t = 3){

  null.stat = c(); null.t = c()

  if(max(order(-stat.list[[1]])[1:3])-min(order(-stat.list[[1]])[1:3]) == 2){
    alt.stat = stat.list[[1]][mean(order(-stat.list[[1]])[1:3])]
    alt.t = mean(order(-stat.list[[1]])[1:3])
  }else if(sum(order(-stat.list[[1]])[1:2] %in% c(1:2))==2){
    alt.stat = stat.list[[1]][1]
    alt.t = 1
  }else if(sum(order(-stat.list[[1]])[1:2] %in% c(length(stat.list[[1]]), (length(stat.list[[1]])-1)))==2){
    alt.stat = stat.list[[1]][length(stat.list[[1]])]
    alt.t = length(stat.list[[1]])
  }else{
    alt.stat = stat.list[[1]][default.t]
    alt.t =  default.t
  }

  for(j in 1:nrow(stat.list[[2]])){
    if(max(order(- stat.list[[2]][j,])[1:3])-min(order(-stat.list[[2]][j,])[1:3]) == 2){
      null.stat[j] = stat.list[[2]][j,mean(order(-stat.list[[2]][j,])[1:3])]
      null.t[j] = mean(order(-stat.list[[2]][j,])[1:3])
    }else if(sum(order(-stat.list[[2]][j,])[1:2] %in% c(1:2))==2){
      null.stat[j] = stat.list[[2]][j,1]
      null.t[j] = 1
    }else if(sum(order(-stat.list[[2]][j,])[1:2] %in% c(length(stat.list[[2]][j,]), (length(stat.list[[2]][j,])-1)) ) == 2){
      null.stat[j] = stat.list[[2]][j,length(length(stat.list[[2]][j,]))]
      null.t[j] = length(length(stat.list[[2]][j,]))
    }else{
      null.stat[j] = stat.list[[2]][j,default.t]
      null.t[j] = default.t
    }
  }

  pval = mean(null.stat >= alt.stat)

  return(list(pvalue = pval, alt.stat = alt.stat, alt.t = alt.t,
              null.stat = null.stat, null.t = null.t))

}

