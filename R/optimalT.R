#' This function empirically searches for data-adaptive optimal scale t for diffusion map embeddings.
#' Input data is a list of statistics under the alternative and those under the null (from permutation).
#' 
#' @param stat.list is interpreted as:
#' \describe{
#'    \item{list of length [2]}{List of a length-[q] vector and [nr x q] matrix comprised of statistics under the alternative and under the null.}
#' }
#' @param default is the position of default Markov time [t] among [q].
#' @return p-value using the opitmal scale; a collection of statistics; optiaml t chosen for the alternative and for [nr] cases of null.
#' @author Youjin Lee
#' @export
#'


print.stat.optimal = function(stat.list, default.t){

  null.stat = c(); null.t = c()
  
  if(max(order(-stat.list[[1]])[1:3])-min(order(-stat.list[[1]])[1:3]) == 2){
    alt.stat = stat.list[[1]][mean(order(-stat.list[[1]])[1:3])]
    alt.t = mean(order(-stat.list[[1]])[1:3])
  }else{
    alt.stat = stat.list[[1]][default.t]
    alt.t = default.t
  }
  
  for(j in 1:nrow(stat.list[[2]])){
    if(max(order(- stat.list[[2]][j,])[1:3])-min(order(-stat.list[[2]][j,])[1:3]) == 2){
      null.stat[j] = stat.list[[2]][j,mean(order(-stat.list[[2]][j,])[1:3])]
      null.t[j] = mean(order(-stat.list[[2]][j,])[1:3])
    }else{
      null.stat[j] = stat.list[[2]][j,default.t]
      null.t[j] = default.t
    }
  }
  
  pval = mean(null.stat >= alt.stat)
  
  return(list(pvalue = pval, alt.stat = alt.stat, alt.t = alt.t, 
              null.stat = null.stat, null.t = null.t))
  
}