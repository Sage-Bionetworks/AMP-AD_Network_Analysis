#' Corr2Cons function 
#' 
#' This function all clustering instances and converts to consensus clustering matrices
#' 
#' @param InstanceList List of all clusterings 


Corr2Cons <- function(InstanceList){
  
  l <- length(InstanceList)
  
  n <- length(InstanceList[[1]])
  
  W <- matrix(0,n,n)
  
  for (i in 1:l){
    
    W <- W + Cluster2Adj(InstanceList[[i]])
    
  }
  
  l2 <- list()
  
  l2$Wp <- W/l 
  l2$Wn <- 1 - W/l
  
  G <- list()
  
  l2$G$V <- c(1:n)
  M <- matrix(1,n,n)
  l2$G$Ep <- (l2$Wp>=l2$Wn)
  l2$G$Em <- (l2$Wp<l2$Wn)
  
  l2$G$Ep[l2$G$Ep] = 1
  l2$G$Em[l2$G$Em] = 1
  
  return(l2)
  
}

#' CC_Pivot
#' 
#' This function performs consensus clustering with CC_Pivot algorithm
#' 
#' @param G A graph object (list) with elements V, Ep and Em 


CC.Pivot <- function(G){
  
  V <- G$V 
  Ep <- G$Ep 
  Em <- G$Em 
  
  if (length(V)<=1){
    return(list(V))
  }
  
  i <- sample(V,1)
  C <- c(i)
  Vp <- c()
  
  
  
  for (j in V){
    if (i == j){
      next
    } 
    
    if (Ep[i,j] == 1){
      C <- c(j,C)
    } else if (Em[i,j] == 1) {
      Vp <- c(j,Vp)
    }
    
  }
  
  #print(C)
  
  n <- sqrt(length(Ep))
  
  G$Ep <- matrix(0,n,n)
  G$Em <- matrix(0,n,n)
  
  G$V <- Vp 
  G$Ep[Vp, Vp] <- Ep[Vp, Vp]
  G$Em[Vp, Vp] <- Em[Vp, Vp]
  
  rm(Ep, Em)
  Gret <- CC.Pivot(G)
  if (is.null(Gret[[1]])){
    return(list(C))
  }
  
  return( c(list(C), Gret))
}  

#' Cluster2Adj 
#' 
#' Converts a cluster to an adjacency matrix 
#' @param ClusterList All clusterings from one clustering instance


Cluster2Adj <- function(ClusterList){
  
  n <- length(ClusterList)
  
  Adj <- outer(1:n, 1:n , FUN=function(r,c) (ClusterList[r] == ClusterList[c])*1)
  # Adj <- matrix(0,n,n)
  # 
  # for (i in 1:n){ 
  #   
  #   for (j in 1:n){ 
  #     
  #     if (ClusterList[i] == ClusterList[j]){
  #       Adj[i,j] <- 1
  #     }
  #     
  #   }
  
  #}
  return(Adj)
}
