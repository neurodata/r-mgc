library(igraph)

notsimple3block = function(popn, p, q, r){
  # generate a 3-block igraph 
  # where within-block edge probabilities are consistent across blocks
  # but between-block edge probabilities are not consistent
  # in a way that most different groups are higher probability
  n.group = 3

  ## 1. generate an outcome variable
  pi = c(1/3, 1/3, 1/3)
  outcome = c()
  for(i in 1:popn){
  	outcome[i] = which(rmultinom(1, 1, pi) == 1)
  }
  
  ## 2. divide sample into three blocks
  ## generate multinomial probability
  multip = matrix(0, nrow = popn, ncol = n.group)
  for(i in 1:popn){
    if(outcome[i] == 1){
      multip[i,] = c(1/2, 1/4, 1/4)
    }else if (outcome[i] == 2){
      multip[i,] = c(1/4, 1/2, 1/4)
    }else{
      multip[i,] = c(1/4, 1/4, 1/2)
    }
  }

  group = rep(0, popn)
  for(i in 1:popn){
    group[i] = which(rmultinom(1, 1, multip[i,]) == 1)
  }
  
  
  ## 3. generate an edge variable
  adj <- matrix(0, popn, popn)
  for(i in 1:popn){
  	for(j in i:popn){
  		if(i == j){
  			adj[i,j] = 0
  		}else{
  			adj[i,j] = rbinom(1, 1, q)
  		}
  		adj[j,i] = adj[i,j]
  	}
  }
  
  prob = matrix(0, popn, popn)
  group1.index = which(group == 1)
  group2.index = which(group == 2)
  group3.index = which(group == 3)
  
  # for group1
  for(i in 1:length(group1.index)){
    for(j in i:length(group1.index)){
      id1 = group1.index[i]
      id2 = group1.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = p
      }  
      adj[id1, id2] = rbinom(1,1,prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]  
    }
  }
  
  # for group2
  for(i in 1:length(group2.index)){
    for(j in i:length(group2.index)){
      id1 = group2.index[i]
      id2 = group2.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = p
      }  
      adj[id1, id2] = rbinom(1,1,prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]  
    }
  }
  
 # for group3
  for(i in 1:length(group3.index)){
    for(j in i:length(group3.index)){
      id1 = group3.index[i]
      id2 = group3.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = p
      }  
      adj[id1, id2] = rbinom(1,1,prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]  
    }
  }
  
    
  # between group1 and group3
  for(i in 1:length(group1.index)){
    for(j in 1:length(group3.index)){
      id1 = group1.index[i]
      id2 = group3.index[j]
      if (id1 == id2){
        prob[id1, id2] = 0.0
      }else{
        prob[id1, id2] = r
      }  
      adj[id1, id2] = rbinom(1,1,prob[id1, id2])
      adj[id2, id1] = adj[id1, id2]  
    }
  }
  
  G = graph.adjacency(adj, "undirected")
 
  V(G)$outcome = outcome
  V(G)$group = group

  return(G)
}  
