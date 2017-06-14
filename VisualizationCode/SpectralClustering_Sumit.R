make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

GetSpectralClusters <- function(A, TopEigen, NoClusters){
  
  k = TopEigen
  k2 = NoClusters 
  
  #generate degree matrix 
  D <- diag(apply(A, 1, sum))
  
  #calculating unnormalized laplacian
  U <- D - A
  
  #finding k smallest eigen vectors
  evL <- eigen(U, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  
  library(stats)
  km <- kmeans(Z, centers=k2, nstart=5)
  
  RetList <- list()
  RetList$km <- km 
  RetList$Z <- Z 
  
  return(RetList)
  
}


