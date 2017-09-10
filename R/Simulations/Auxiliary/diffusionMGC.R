NetworkTest.diffusion.stat = function(G, X, option, diffusion, t.range, n.perm, q){
  
  if(!is.connected(G)){
    G <- connect_graph(G)
  }
  
  A <- as.matrix(get.adjacency(G))
  P <- A / rowSums(A)
  
  if(diffusion == FALSE){ 
    Dx = as.matrix(dist((A)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      result = MGCPermutationTest(Dx, Dy, n.perm, option = 'mgc')
      return(result)
    }else if(option == 2){
      pvalues = dcov.test(Dx, Dy, index = 1.0, R = n.perm)$p.value
      return(pvalues)
    }else if(option == 3){
      pvalues = hhg.test(Dx, Dy, nr.perm = n.perm)$perm.pval.hhg.sl
      return(pvalues)
    }
  }

  U.list <- dmap.q(P, t.range, q) 
  
  mgc.stat = c(); dcor.stat = c(); hhg.stat = c()
  
  for(s in 1:length(t.range)){
    
    U <- U.list[[s]]
    Dx = as.matrix(dist((U)), diag = TRUE, upper = TRUE) 
    Dy = as.matrix(dist((X)), diag = TRUE, upper = TRUE)
    
    if(option == 1){
      mgc = MGCSampleStat(Dx, Dy, option = 'mgc')
      mgc.stat[s] = mgc$statMGC
    }
    
    if(option == 2){
      dcor = dcor.ttest(Dx, Dy, distance = TRUE)
      dcor.stat[s] = dcor$statistic
    }
    
    
    if(option == 3){
      hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
      hhg.stat[s] = hhg$sum.lr
    }
  }

  mgc.pval = 0; dcor.pval = 0; hhg.pval = 0; tmp = matrix(0, nrow = n.perm, ncol = length(t.range))

  for(r in 1:n.perm){

    for(s in 1:length(t.range)){
      U <- U.list[[s]]
      per=sample(length(X));
      newX = X[per]
      Dx = as.matrix(dist(U), diag = TRUE, upper = TRUE) 
      Dy = as.matrix(dist(newX), diag = TRUE, upper = TRUE)
    
      if(option == 1){
        mgc = MGCSampleStat(Dx, Dy, option = 'mgc')
        tmp[r,s] = mgc$statMGC
      }
    
      if(option == 2){
        dcor = dcor.ttest(Dx, Dy, distance = TRUE)
        tmp[r,s] = dcor$statistic
      }
    
      if(option == 3){
        hhg = hhg.test(Dx, Dy, nr.perm = n.perm)
        tmp[r,s] = hhg$sum.lr
      }
    } 
  }

  if(option == 1) return(list(mgc.stat, tmp))

  if(option == 2) return(list(dcor.stat, tmp))

  if(option == 3) return(list(hhg.stat, tmp))

 
}


