library(DTMCPack)
library(Matrix)

#### print out diffusion maps from undirected grph
dmap = function(trans, times){
    # require : DTMCPack, Matrix
    # trans : a transition probability
    # times : a set of diffusion times you want to print out
    maps = list()
    
    # stationary distribution (probability)
    pi = statdistr(trans)
    pi.mat1 = as.matrix(Diagonal(length(pi), pi))
    pi.mat2 = as.matrix(Diagonal(length(pi), pi^(-1/2)))
    
    # symmetric kernel
    Q = sqrt(pi.mat1) %*% trans %*% pi.mat2
    
    lambda = Re(eigen(Q)$values)[Im(eigen(Q)$value) == 0]
    psi = Re(eigen(Q)$vectors)[,Im(eigen(Q)$value) == 0]
    
    phi = pi.mat2 %*% psi
    
    Lambda = as.matrix(Diagonal(length(lambda), lambda))
    
    # each row of maps is a diffusion maps of each vertex
    for(t in 1:length(times)){
        maps[[t]] = phi %*% Lambda^(times[t])  
    }

    return(maps)
}
