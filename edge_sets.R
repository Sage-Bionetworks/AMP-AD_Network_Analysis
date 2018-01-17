#pull edge lists
synapser::synLogin()

#pull modules from synapse
aggregateModules <- synapser::synTableQuery("select * from syn11182793")
aggMods <-aggregateModules$asDataFrame()

#pull network from synapse
foo <- synapser::synTableQuery('select * from syn8681664 where method=\'bic\' and assay =\'RNAseq\'')
foo <- foo$asDataFrame()

loadBic <- function(synId){
  foo<-synapser::synGet(synId)
  load(foo$path)
  net <- as.matrix(bicNetworks$network)
  net <- net+t(net)
  return(net)
}

ampNetworks <- lapply(foo$id,loadBic)
names(ampNetworks) <- foo$tissueTypeAbrv

#convert to edgelists
foobar <- lapply(ampNetworks,igraph::graph_from_adjacency_matrix,mode='undirected')
foobar2 <- lapply(foobar,igraph::as_data_frame)

dummyFun <- function(x,y){
  genes <- unique(c(x$from,x$to))
  mapTable <- NA
  while(is.na(mapTable)){
    try(mapTable <- utilityFunctions::convertEnsemblToHgnc(genes),silent=T)
  }
  x <- dplyr::left_join(x,mapTable,by=c('from'='ensembl_gene_id'))
  x<- dplyr::left_join(x,mapTable,by=c('to'='ensembl_gene_id'))
  colnames(x) <- c('geneA_ensembl_gene_id',
                   'geneB_ensembl_gene_id',
                   'geneA_external_gene_name',
                   'geneB_external_gene_name')
  x$brainRegion <- rep(y,nrow(x))
  return(x)
}

foobar3 <- mapply(dummyFun,
                  foobar2,
                  names(foobar2),
                  SIMPLIFY=FALSE)
edgeListManifest <- do.call(rbind,foobar3)
#foobar3 <- 
#push to synapse with synapser

permLink <- githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'module_ranking',
                               repositoryPath = 'edge_sets.R')

edgeListObj<-rSynapseUtilities::pushDf2Synapse(df = edgeListManifest,
                                  fileName = 'metanetwork_edges_amp_ad.csv',
                                  synapseFolderId = 'syn7525089',
                                  annos = c('analysisType' = 'statisticalNetworkReconstruction',
                                            'assay' = 'RNAseq',
                                            'consortium' = 'AMP-AD',
                                            'dataSubType' = 'geneExp',
                                            'dataType' = 'analysis',
                                            'fileType' = 'csv',
                                            'method' = 'bic',
                                            'normalizationStatus' = 'TRUE',
                                            'normalizationType' = 'CPM',
                                            'organism' = 'HomoSapiens',
                                            'summaryLevel' = 'gene'),
                                  comment = 'first version of edge list for wall of targets visulization',
                                  usedVector = foo$id,
                                  executedVector = permLink,
                                  activityName1 = 'Make Edge List',
                                  activityDescription1 = 'pull bic networks, turn into edge list, and add human readable names')

