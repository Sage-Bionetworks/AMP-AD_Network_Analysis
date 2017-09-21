computeBrainRegionConsensus <- function(brainRegion,
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
  #synapseClient::synapseLogin()
  
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
  
  #print(head(partition.adj))
  # Randomise gene order
  partition.adj <- partition.adj[sample(1:dim(partition.adj)[1], dim(partition.adj)[1]), ]
  mod <- metanetwork::findModules.consensusCluster(d = t(partition.adj), 
                                                   maxK = 100, 
                                                   reps = 50, 
                                                   pItem = 0.8, 
                                                   pFeature = 1,
                                                   clusterAlg = "hc", 
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
  parallel::stopCluster(cl)
  return(mod)
}

pushToSynapseWrapper <- function(df,
                                 fileName,
                                 synapseFolderId,
                                 annos,
                                 comment,
                                 usedVector,
                                 executedVector,
                                 activityName1,
                                 activityDescription1){
  synapseClient::synapseLogin()
  
  #write df to file
  write.csv(df,file=fileName,quote=F)
  
  #make File object (where it goes on synapse and comment)
  foo <- synapseClient::File(fileName,
                             parentId=synapseFolderId,
                             versionComment=comment)
  
  #apply annotations
  synapseClient::synSetAnnotations(foo) = as.list(annos)
  
  #push to synapse
  foo = synapseClient::synStore(foo,
                                used = as.list(usedVector),
                                executed = as.list(executedVector),
                                activityName = activityName1,
                                activityDescription = activityDescription1)
  
  #return the synapse object
  return(foo)
}


ClusterLowFreqConsensus <- function(Labs,cutoff){
  
  library(plyr)
  m <- max(Labs) + 1 
  Lab_In <- which(plyr::count(Labs)$freq <= cutoff)
  In <- which(Labs %in% Lab_In)
  Labs[In] <- m + 1 
  return(Labs)
}


GenConsMods <- function(Module,mthd,ModuleList,br){
  
  l <- list()
  external_gene_name <- unique(ModuleList$external_gene_name)
  brainRegion <- c(rep(br,length(external_gene_name)))
  method <- c(rep(mthd,length(external_gene_name)))
  GeneID <- c()
  
  for (i in 1:length(external_gene_name)){
    In <- which(ModuleList$external_gene_name %in% external_gene_name[i])
    GeneID <- c(GeneID, ModuleList$GeneID[In[1]])
  }
  
  ModuleName <- c(rep(br,length(external_gene_name)))
  ModuleNameFull <- c(rep(br,length(external_gene_name)))
  
  for (i in 1:length(external_gene_name)){
    ModuleName[i] <- paste0(mthd,Module[i])
    ModuleNameFull[i] <- paste0(ModuleName[i],br)
  }
  
  #l <- list(GeneID, Module, method, ModuleName, external_gene_name,
  #          brainRegion, ModuleNameFull)
  l$GeneID <- GeneID
  l$Module <- Module
  l$method <- method
  l$ModuleName <- ModuleName
  l$external_gene_name <- external_gene_name
  l$brainRegion <- brainRegion
  l$ModuleNameFull <- ModuleNameFull
  
  return(data.frame(l, stringsAsFactors = F))
}