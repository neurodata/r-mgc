library(igraph)
source("j-RDPG.R")

n <- 50
m <- 10
d <- 2

Z = matrix(0, ncol = 2, nrow = n)
Z[1:(n/2),1] = 1
Z[(n/2+1):n,2] = 1

V = svd(Z)$u

sample_graph <- function() {
  p1 = runif(1)
  p2 = runif(1)
  q = runif(1)
  pm <- cbind( c(p1, q), c(q, p2))
  g <-sample_sbm(n, pref.matrix=pm, block.sizes=c(n/2,n/2))
  A = get.adjacency(g)
  A
}

Adj_list <- sapply(1:m, function(i) sample_graph())

jrdpg <- joint_latent_positions(Adj_list = Adj_list, d = 2,K = d, ASE = FALSE, 
                                par = FALSE)


plot(jrdpg$V)
subspace_distance(V, jrdpg$V)

