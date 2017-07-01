computeBrainRegionConsensus(brainRegion,
                            manifestId,
                            seed=5){
  set.seed(seed)
  
  nc = parallel::detectCores()
  if (nc > 2){
    cl = parallel::makeCluster(nc - 2)
  } else {
    cl = parallel::makeCluster(1)
  }
  doParallel::registerDoParallel(cl)
  
  #### Login to synapse ####
  synapseClient::synapseLogin()
  
  #### Get individual partitions from synapse and formulate partition adjacency matrix ####
  queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",brainRegion,"' ) )")
  indMods <- synapseClient::synTableQuery(queryString)@values
  
  convertFormat <- function(methodVal,modTab){
    foo <- dplyr::filter(modTab,method==methodVal)
    bar <- dplyr::select(foo,GeneID,Module)
    colnames(bar) <- c('Gene.ID','moduleNumber')
    return(bar)
  }
  
  partition.adj <- lapply(unique(indMods$method),convertFormat,indMods)
  names(partition.adj) <- unique(indMods$method)
  library(dplyr)
  partition.adj <- mapply(function(mod, method){
    mod = mod %>%
      dplyr::select(Gene.ID, moduleNumber) %>%
      dplyr::mutate(value = 1,
                    moduleNumber = paste0(method,'.',moduleNumber)) %>%
      tidyr::spread(moduleNumber, value)
  }, partition.adj, names(partition.adj), SIMPLIFY = F) %>%
    plyr::join_all(type="full")
  partition.adj[is.na(partition.adj)] <- 0
  rownames(partition.adj) <- partition.adj$Gene.ID
  partition.adj$Gene.ID <- NULL
  
  # Randomise gene order
  partition.adj <- partition.adj[sample(1:dim(partition.adj)[1], dim(partition.adj)[1]), ]
  mod <- metanetwork::findModules.consensusCluster(d = t(partition.adj[1:400,]), 
                                                   maxK = 100, 
                                                   reps = 50, 
                                                   pItem = 0.8, 
                                                   pFeature = 1,
                                                   clusterAlg = "kmeans", 
                                                   innerLinkage = "average", 
                                                   distance = "pearson",
                                                   changeCDFArea = 0.001, 
                                                   nbreaks = 10, 
                                                   seed = 1,
                                                   weightsItem = NULL, 
                                                   weightsFeature = NULL, 
                                                   corUse = "everything",
                                                   verbose = F)
  
  # Find modularity quality metrics
  mod <- data.frame(mod,stringsAsFactors=F)
  stopCluster(cl)
  return(mod)
}
dlpfcConsensus <- computeBrainRegionConsensus('DLPFC','syn10146524')
tcxConsensus <- computeBrainRegionConsensus('TCX','syn10146524')
cbeConsensus <- computeBrainRegionConsensus('CBE','syn10146524')
fpConsensus <- computeBrainRegionConsensus('FP','syn10146524')
ifgConsensus <- computeBrainRegionConsensus('IFG','syn10146524')
phgConsensus <- computeBrainRegionConsensus('PHG','syn10146524')
stgConsensus <- computeBrainRegionConsensus('STG','syn10146524')