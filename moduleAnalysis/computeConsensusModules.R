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
                                                   verbose = F,
                                                   useParallelFlag = T)
  
  # Find modularity quality metrics
  mod <- data.frame(mod,stringsAsFactors=F)
  parallel::stopCluster(cl)
  return(mod)
}
dlpfcConsensus <- computeBrainRegionConsensus('DLPFC','syn10146524')
write.csv(dlpfcConsensus,file='dlpfc.csv',quote=F)
tcxConsensus <- computeBrainRegionConsensus('TCX','syn10146524')
write.csv(tcxConsensus,file='tcx.csv',quote=F)
cbeConsensus <- computeBrainRegionConsensus('CBE','syn10146524')
write.csv(cbeConsensus,file='cbe.csv',quote=F)
fpConsensus <- computeBrainRegionConsensus('FP','syn10146524')
write.csv(fpConsensus,file='fp.csv',quote=F)
ifgConsensus <- computeBrainRegionConsensus('IFG','syn10146524')
write.csv(ifgConsensus,file='ifg.csv',quote=F)
phgConsensus <- computeBrainRegionConsensus('PHG','syn10146524')
write.csv(phgConsensus,file='phg.csv',quote=F)
stgConsensus <- computeBrainRegionConsensus('STG','syn10146524')
write.csv(stgConsensus,file='stg.csv',quote=F)

#upload to synapse
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
###read in from files into a list
fileNames <- c('dlpfc.csv',
               'cbe.csv',
               'tcx.csv',
               'fp.csv',
               'ifg.csv',
               'phg.csv',
               'stg.csv')

#change working directory
setwd('~/AMP-AD_Network_Analysis/moduleAnalysis/')
###modify so that it looks right
modList <- lapply(fileNames,
                  data.table::fread,
                  data.table=F)
modList <- lapply(modList,function(x) dplyr::select(x,-V1))
###push to synapse using the above function
mapply(pushToSynapseWrapper,
       modList,
       c('consensusDLPFC.csv',
         'consensusCBE.csv',
         'consensusTCX.csv',
         'consensusFP.csv',
         'consensusIFG.csv',
         'consensusPHG.csv',
         'consensusSTG.csv'),
       as.list(rep('syn10158370',7)),
       list(c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='dorsolateralPrefrontalCortex',study='ROSMAP',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='DLPFC',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='cerebellum',study='MayoRNAseq',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='CBE',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='temporalCortex',study='MayoRNAseq',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='TCX',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='frontalPole',study='MSBB',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='FP',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='inferiorFrontalGyrus',study='MSBB',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='IFG',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='parahippocampalGyrus',study='MSBB',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='PHG',geneIdentifiers='Ensembl',method='consensusKmeans'),
            c('analysisType'='moduleIdentification',fileType='csv',assay='RNAseq',tissueOfOrigin='superiorTemporalGyrus',study='MSBB',summaryLevel='gene',organism='HomoSapiens',consortium='AMP-AD',dataType='analysis',tissueTypeAbrv='STG',geneIdentifiers='Ensembl',method='consensusKmeans')),
       list('consensus module generation for rosmap dlpfc',
            'consensus module generation for mayornaseq cbe',
            'consensus module generation for mayornaseq tcx',
            'consensus module generation for msbb fp',
            'consensus module generation for msbb ifg',
            'consensus module generation for msbb phg',
            'consensus module generation for msbb stg'),
       as.list(rep('syn10146524',7)),
       rep(list(c('https://github.com/Sage-Bionetworks/AMP-AD_Network_Analysis/blob/dd8a114ea8e60b24efa17b5017724caabf07edb3/moduleAnalysis/computeConsensusModules.R',
                  'https://github.com/Sage-Bionetworks/metanetwork/blob/0e7a52e9401c9979632faf475fb3d9ad0249736c/R/findModules.consensusCluster.R')),7),
       as.list(rep('consensus module construction',7)),
       as.list(rep('consensus module building workflow',7)))
###add a column for consensus
dfTest <- mapply(function(df,
                          brainRegion){
  foo <- dplyr::select(df,Gene.ID,moduleLabel)
  colnames(foo) <- c('GeneID','Module')
  foo$method <- rep('consensus',nrow(foo))
  foo$ModuleName <- paste0(foo$method,foo$Module)
  bar <- utilityFunctions::convertEnsemblToHgnc(foo$GeneID)
  foo <- dplyr::left_join(foo,bar,by=c('GeneID'='ensembl_gene_id'))
  foo$brainRegion <- rep(brainRegion,nrow(foo))
  foo$ModuleNameFull <- paste0(foo$ModuleName,foo$brainRegion)
  return(foo)
},
modList,
as.list(c('DLPFC',
          'CBE',
          'TCX',
          'FP',
          'IFG',
          'PHG',
          'STG')),
SIMPLIFY=FALSE)
###make schema correct, combine into a single data frame
fullConsensus <- do.call(rbind,dfTest)
rSynapseUtilities::makeTable(fullConsensus,'consensus module manifest July 13 2017','syn2370594')

###upload to a consensus table for downstream analyses

