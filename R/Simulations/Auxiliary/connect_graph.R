library(igraph)

connect_graph = function(G){
  # this function turns disconnected graph into connected graph
  # by finding components, chosing one vertex for each 
  # and create an edge between them
  
  connected_G = G
  V(connected_G)$name = c(1:length(V(G)))
  
  if (!is.connected(connected_G)){ # if there exists more than one component
    
    num = no.clusters(connected_G) # number of clusters
    
    selected_node = c() # the index vector selected from each cluster
    
    for (i in 1:num){
      
      sub = induced.subgraph(connected_G, clusters(connected_G)$membership == i)
      
      if(length(V(sub)$name) == 1){
        selected_node[i] = V(sub)$name
      }else{
        selected_node[i] = sample(V(sub)$name, 1) # select one sample from sub
      }
    }
    
    for (j in 2:num){
      connected_G = connected_G + 
        edges(c(V(connected_G)[V(connected_G)$name == selected_node[1]] ,
                V(connected_G)[V(connected_G)$name == selected_node[j]]))
      }
  }
  
  if(is.connected(connected_G)){
    return(connected_G)
  }else{
    return(NULL)
  }
  
}
