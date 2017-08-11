ConvertModsToMat <- function(ModuleList){ 
  
  GeneList <- unique(ModuleList$external_gene_name)
  MethodList <- unique(ModuleList$method)
  N <- length(GeneList)
  
  InstanceList <- list()
  #create cluster label representation
  for (i in 1:length(MethodList)){
    Clust <- rep(0,N)
    
    In0 <- which(ModuleList$method %in% MethodList[i])
    UnqMods <- unique(ModuleList$ModuleName[In0])
    for( j in 1:length(UnqMods)){
      
      In <- which(ModuleList$ModuleName %in% UnqMods[j])
      #print(head(In))
      In2 <- which(GeneList %in% ModuleList$external_gene_name[In])
      #print(head(In2))
      Clust[In2] <- j 
      
    }
    InstanceList[[i]] <- Clust
  }
  
  return(InstanceList)
  
  }