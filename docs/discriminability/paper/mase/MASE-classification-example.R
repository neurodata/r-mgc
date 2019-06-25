source("R/loadAll.R")
# Joint embedding of graphs code
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/jesusdaniel/JEG/master/Code/joint_embedding.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
# Multiple RDPG
library(multiRDPG)
# Omnibus embedding
source("R/omnibus-embedding.R")

### Simulate a sample of four classes of graphs
m = 40
n = 200
d=2
# Separation between classes
alpha = 1


# plot eigenvectors on 2 dimensions
plot_eig <- function(V) {
  cols <- c(rep("orange", n/2),rep("red", n/2))
  xlims = c(-0.2, 0.2)
  ylims = c(-0.2, 0.2)
  df <- data.frame(x = V[,1], y = V[,2], col=cols)
  print(ggplot(df, aes(x=x, y=y)) + 
          #geom_hline(yintercept = 0)+
          #geom_vline(xintercept = 0)+
          geom_point(aes(col="black",fill=cols),shape=21,size=3) +
          xlim(xlims) + ylim(ylims) + xlab("") +
          ylab("") +
          scale_color_manual(values=c("black"))+
          scale_fill_manual(values=c("orange", "red"))+
          theme_bw() + theme(legend.position="none")) 
}




# Simulate four classes of B matrices
Bs = lapply(1:m, function(i) {
  B0 = matrix(0.15,  d, d)
  b1 <- if(i%%4==1) {
    c(0.1, 0.0, 0.1)
  } else{ if(i%%4==2) {
    -c(0.1, 0, 0.1)
  } else{ if(i%%4==3) {
    c(0.1, 0.0, 0.0)
  }else {
    c(0.0, 0.0, 0.1)
  }}}
  B = matrix(c(b1[1:2], b1[2], b1[3]), nrow = 2)
  B0 + alpha * B
})
# Class labels
labels <- factor(rep(1:4, m/4))

# Memberships of a two-block SBM
Z <- matrix(0, n, 2)
Z[1:(n/2), 1] = Z[(n/2 + 1):n,2] = 1

# Generate sample
Adj_list <- lapply(Bs, function(B) {
  get.adjacency(sample_sbm(n, B, block.sizes = rep(n/2, 2),
                           directed = FALSE))
})
# Convert to non-sparse matrices
Adj_list_m <- lapply(Adj_list, as.matrix)

plot_adjmatrix(Adj_list[[1]])
plot_adjmatrix(Adj_list[[2]])
plot_adjmatrix(Adj_list[[3]])
plot_adjmatrix(Adj_list[[4]])


# MASE
mase.embed <- mase(Adj_list, d = 2)
# 2D plot of latent positions
plot_eig(jrdpg$V)


# JEG - random initialization
jeg <- multidembed(Adj_list, d = 3, Innitialize = 2)

# MRDPG
mrdpg <- multiRDPG(Adj_list_m, d = 3)

# OMNI embedding
omni <- OMNI_matrix(Adj_list, 2)




# Obtain graph coordinates
mase.R.vector = t(sapply(mase.embed$R, function(R) R[upper.tri(R, diag=TRUE)]))

# MRDPG coordinates
mrdpg_L = t(sapply(mrdpg$Lambda, diag))



# 3d plots of the network coordinates
plot.3d((mase.R.vector), col = labels)
plot.3d((jeg$lambda), col = labels)
plot.3d((mrdpg_L), col = labels)

# Observe that MASE splits the classes accurately!



##classification accuracy

# create test set
Adj_list_test <- lapply(Bs, function(B) {
  get.adjacency(sample_sbm(n, B, block.sizes = rep(n/2, 2),
                           directed = FALSE))
})
labels_test = labels

Bhat_test <- project_networks(Adj_list_test, mase.embed$V)

# use knn
Xtrain <- t(sapply(mase.embed$R, as.vector))
Xtest <- t(sapply(Bhat_test, as.vector))

library(class)
# MASE classification accuracy
sum(knn(train = Xtrain, test = Xtest, cl = labels, k = 1) == labels) /  length(labels)


# mrdpg
# project mrdpg test set
Xtrain.mrdpg <- mrdpg_L
Xtest.mrdpg <- t(sapply(Adj_list_test, function(A) diag((t(mrdpg$U) %*% crossprod(A, mrdpg$U) ))))
sum(knn(train = Xtrain.mrdpg, test = Xtest.mrdpg, cl = labels, k = 1) == labels) /  length(labels)



# jeg
Xtrain.jeg <- jeg$lambda
Xtest.jeg <- multiembed_test(A = Adj_list_test, h = jeg$h)
sum(knn(train = Xtrain.jeg, test = Xtest.jeg, cl = labels, k = 1) == labels) /  length(labels)


