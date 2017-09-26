#' This function tests the equivalence of K distributions using distance-based test including MGC.
#' Input data from K populations are first transformed into pooled data Z and label vector Y.
#' then implements test statistics upon distance matrices from each data X and Y.
#'
#' @param dat is interpreted as:
#' \describe{
#'    \item{list of length [k]}{dat is a k-length list of the same size of vectors or matrices}
#' }
#' @param dist='Euclidean' is a string that specifies distance matrix to apply.
#' @return result : list of results of mgc.stat testing equivalence of K-sample test.
#' @author Youjin Lee
#' @export
#'

mgc.ksample = function(dat, dist = 'Euclidean'){
    if(class(dat) == 'list'){
        K = length(dat)
    }else{
        print("Data set should be in the list format")
    }

    if(class(dat[[1]]) == 'numeric'){
        X = c(); Y = c()
        for(i in 1:K){
            X = c(X, dat[[i]])
            Y = c(Y, rep(i, length(dat[[i]])))
        }
        Dx = as.matrix(dist(X, diag = TRUE, upper = TRUE))
        Dy = as.matrix(dist(Y, diag = TRUE, upper = TRUE))
        Dy = ifelse(Dy > 0, 1, 0)
    }

    if(class(dat[[1]]) == 'matrix'){
         X = c(); Y = c()
        for(i in 1:K){
            X = rbind(X, dat[[i]])
            Y = c(Y, rep(i, nrow(dat[[i]])))
        }
        Dx = as.matrix(dist(X, diag = TRUE, upper = TRUE))
        Dy = as.matrix(dist(Y, diag = TRUE, upper = TRUE))
    }

    result =  mgc.test(Dx, Dy, rep=1000, option = 'mgc')

    return(result)
}
