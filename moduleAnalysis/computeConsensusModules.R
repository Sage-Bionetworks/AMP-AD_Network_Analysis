exprId <- 'syn8303260'
brainRegion <- 'DLPFC'
manifestId <- 'syn10146524'

computeBrainRegionConsensus(exprId,brainRegion,manifestId,seed=5){
  set.seed(seed)
  #### Load Libraries ####
  # library(data.table)
  # library(tidyr)
  # library(plyr)
  # library(dplyr)
  # 
  # library(igraph)
  # library(metanetwork)
  # 
  # library(synapseClient)
  # library(githubr)
  # 
  # library(parallel)
  # library(doParallel)
  # library(foreach)
  
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
  partition.adj = mapply(function(mod, method){
    mod = mod %>%
      dplyr::select(Gene.ID, moduleNumber) %>%
      dplyr::mutate(value = 1,
                    moduleNumber = paste0(method,'.',moduleNumber)) %>%
      tidyr::spread(moduleNumber, value)
  }, partition.adj, names(partition.adj), SIMPLIFY = F) %>%
    plyr::join_all(type="full")
  partition.adj[is.na(partition.adj)] = 0
  rownames(partition.adj) = partition.adj$Gene.ID
  partition.adj$Gene.ID = NULL
  
  # Randomise gene order
  partition.adj = partition.adj[sample(1:dim(partition.adj)[1], dim(partition.adj)[1]), ]
  
  #### Compute consensus modules using specified algorithm ####
  # Compute consensus modules nreps and choose the best solution
  
  #library('ConsensusClusterPlus')
  a1 <- ConsensusClusterPlus::ConsensusClusterPlus(d = t(partition.adj[1:300,]), 
                                             maxK = 10,
                                             reps = 10,
                                             pItem = 0.8,
                                             pFeature = 1,
                                             clusterAlg = "km",
                                             innerLinkage = "average",
                                             distance = "euclidean",
                                             seed = 1,
                                             weightsItem = NULL,
                                             weightsFeature = NULL,
                                             corUse = "everything",
                                             verbose = F)
  
  
  mod <- metanetwork::findModules.consensusCluster(d = t(partition.adj), maxK = 100, reps = 1, pItem = 0.8, pFeature = 1,
                                                   clusterAlg = "kmeans", innerLinkage = "average", distance = "pearson",
                                                   changeCDFArea = 0.001, nbreaks = 10, seed = 1,
                                                   weightsItem = NULL, weightsFeature = NULL, corUse = "everything",
                                                   verbose = F)
  
  # Find modularity quality metrics
  mod = as.data.frame(mod)
  rownames(mod) = mod$Gene.ID
  mod = mod[rownames(adj),]
  NQ = metanetwork::compute.LocalModularity(adj, mod)
  Q = metanetwork::compute.Modularity(adj, mod, method = 'Newman1')
  Qds = metanetwork::compute.ModularityDensity(adj, mod)
  module.qc.metrics = metanetwork::compute.ModuleQualityMetric(adj, mod)
  
  #### Store results in synapse ####
  # Create a modules folder
  fold = Folder(name = 'Modules', parentId = synGet(bic.obj@properties$parentId, downloadFile = F)@properties$parentId)
  fold = synapseClient::synStore(fold)
  
  # Create a consensus modules folder
  fold1 = Folder(name = 'consensus_kmeans', parentId = fold@properties$id)
  fold1 = synStore(fold1)
  
  # Write results to synapse
  system(paste('mkdir',bicNet.id))
  write.table(mod, file = paste0(bicNet.id,'/Consensus.',cons.method,'.',run.id,'.modules.tsv'), row.names=F, quote=F, sep = '\t')
  obj = synapseClient::File(paste0(bicNet.id,'/Consensus.',cons.method,'.',run.id,'.modules.tsv'), parentId = fold1$properties$id)
  synapseClient::annotations(obj) = synapseClient::annotations(bic.obj)
  obj$annotations$fileType = "tsv"
  obj$annotations$analysisType = "consensusModuleIdentification"
  obj$annotations$method = cons.method
  obj$annotations$Q = Q
  obj$annotations$NQ = NQ
  obj$annotations$Qds = Qds
  obj$annotations$columnScaled = TRUE
  obj$annotations$winsorized = TRUE 
  obj$annotations$deprecated = FALSE
  obj = synapseClient::synStore(obj, used = all.used.ids, executed = thisFile, activityName = 'Consensus Module Identification')
  
  write.table(module.qc.metrics, file = paste0(bicNet.id,'/Consensus.',cons.method,'.',run.id,'.moduleQCMetrics.tsv'), row.names=F, quote=F, sep = '\t')
  obj.qc = synapseClient::File(paste0(bicNet.id,'/Consensus.',cons.method,'.',run.id,'.moduleQCMetrics.tsv'), parentId = fold1$properties$id)
  synapseClient::annotations(obj.qc) = synapseClient::annotations(obj)
  obj.qc$annotations$analysisType = "moduleQC"
  obj.qc = synapseClient::synStore(obj.qc, activity = synGetActivity(obj))
  
  stopCluster(cl)

}