####login to Synapse
synapseClient::synapseLogin()

####pull top rankings from table
buildTargetedModules <- function(tissueType){
  library(dplyr)
  bar <- synapseClient::synTableQuery("SELECT * FROM syn10879238")@values %>%
    dplyr::filter(ModuleMethod!='consensus' & ModuleBrainRegion==tissueType)
  allMods <- synapseClient::synTableQuery("SELECT * FROM syn10338156")@values
  pairwiseString <- paste0("SELECT * FROM syn10339153 where ModuleNameFull like \'%",
                           tissueType,
                           "\' and category like \'%",
                           tissueType,
                           "\'")
  pairwise <- synapseClient::synTableQuery(pairwiseString)@values %>%
    utilityFunctions::removeSwappedDupKeyValueDf() %>% 
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%bar$ModuleNameFull & to %in% bar$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)
  
  res <- list()
  res$moduleGraph <- pairwise
  
  graph1 <- igraph::graph_from_data_frame(res$moduleGraph,directed=FALSE)
  test1 <- igraph::optimal.community(graph1)
  metaGraph <- data.frame(ModuleNameFull = test1$names,
                          metaModule = test1$membership,
                          stringsAsFactors=F)

  res$moduleGraphCommunities <- metaGraph
  
  getMajority <- function(df){
    masterTableHGNC <- table(df$external_gene_name)
    masterTableENSG <- table(df$GeneID)
    #masterTable <- masterTable/len1
    masterTableHGNC <- which(masterTableHGNC > 1)
    masterTableENSG <- which(masterTableENSG > 1)
    gen <- list()
    gen$ensg <- names(masterTableENSG)
    gen$hgnc <- names(masterTableHGNC)
    return(gen)
  }
  
  mods <- unique(res$moduleGraphCommunities$metaModule)
  dfList <- lapply(mods,function(x,df,allMods){
    library(dplyr)
    dplyr::filter(allMods,ModuleNameFull %in% df$ModuleNameFull[df$metaModule==x]) %>%
      return},
    res$moduleGraphCommunities,
    allMods)

  res$mods <- lapply(dfList,getMajority)
  names(res$mods) <- paste0(tissueType,WGCNA::labels2colors(mods))
  list2Df <- function(x,module,tissueType){
    mod <- data.frame(GeneID = x$ensg,
                      Module = rep(module,length(x$ensg)),
                      method = rep('aggregate',length(x$ensg)),
                      ModuleName = paste0('aggregate',rep(module,length(x$ensg))),
                      brainRegion = rep(tissueType,length(x$ensg)),
                      ModuleNameFull = paste0('aggregate',rep(module,length(x$ensg)),tissueType),
                      stringsAsFactors=F)
    return(mod)}
  res$df<-mapply(list2Df,
                 x = res$mods,
                 module = names(res$mods),
                 MoreArgs = list(tissueType=tissueType),
                 SIMPLIFY=F)
  res$df <- do.call(rbind,res$df)
  exg <- utilityFunctions::convertEnsemblToHgnc(res$df$GeneID)
  res$df <- dplyr::left_join(res$df,exg,by = c('GeneID' = 'ensembl_gene_id'))
  return(res)
}
DLPFCres <- buildTargetedModules('DLPFC')
CBEres <- buildTargetedModules('CBE')
TCXres <- buildTargetedModules('TCX')
IFGres <- buildTargetedModules('IFG')
STGres <- buildTargetedModules('STG')
PHGres <- buildTargetedModules('PHG')
FPres <- buildTargetedModules('FP')

AggregateModuleManifest <- rbind(DLPFCres$df,
                                 CBEres$df,
                                 TCXres$df,
                                 IFGres$df,
                                 STGres$df,
                                 PHGres$df,
                                 FPres$df)

rSynapseUtilities::makeTable(AggregateModuleManifest,'collapsed ad meta modules september 21 2017',projectId='syn5569099')


###combine modules into a single list

foo <- c(DLPFCres$mods,
         CBEres$mods,
         TCXres$mods,
         IFGres$mods,
         STGres$mods,
         PHGres$mods,
         FPres$mods)

###
foo <- lapply(foo,function(x){return(x$hgnc)})

genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))


system.time(aaaa <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                                  foo,
                                                  GeneSets$Cell_Markers,
                                                  unique(unlist(foo))))

####pull expression data-set
source('dataPulling/pullExpressionAndPhenoWinsorized.R')

####run consensus clustering per tissue type via interesting modules
tissue <- 'DLPFC'

####pull pairwise relationships between modules
pairwise <- synapseClient::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like \'%DLPFC\' and category like \'%DLPFC\'")@values
pairwise <- utilityFunctions::removeSwappedDupKeyValueDf(pairwise)
pairwise <- dplyr::mutate(pairwise,adj=p.adjust(fisherPval,method='bonferroni'))
pairwise <- dplyr::filter(pairwise,adj<=0.05)


keep_mods <- dplyr::filter(bar,ModuleBrainRegion==tissue)
View(keep_mods)
pairwise <- dplyr::filter(pairwise,from%in%keep_mods$ModuleNameFull & to %in% keep_mods$ModuleNameFull)





pairwise <- dplyr::mutate(pairwise, weight = 1/fisherOR)
write.csv(pairwise,file='~/Desktop/dlpfc_targeted.csv',quote=F)


graph1 <- igraph::graph_from_data_frame(pairwise,directed=FALSE)

test1 <- igraph::optimal.community(graph1)

metaGraph <- data.frame(ModuleNameFull = test1$names,
                        metaModule = test1$membership,
                        stringsAsFactors=F)
write.csv(metaGraph,file='~/Desktop/metaGraphDLPFC.csv',quote=F)

modTargeted1 <- dplyr::filter(allMods,ModuleNameFull %in% metaGraph$ModuleNameFull[metaGraph$metaModule==1])

###voting - module built based on being seen in majority of modules in super cluster
getMajority <- function(df){
  masterTable <- table(df$external_gene_name)
  len1 <- length(unique(df$ModuleNameFull))
  #masterTable <- masterTable/len1
  masterTable <- which(masterTable > 1)
  return(names(masterTable))
}


keep1a <- getMajority(modTargeted1)
pheatmap::pheatmap(cor(geneExpressionForAnalysis$rosmapDLPFC[,keep1a]))

View(modTargeted)

library(dplyr)
genes <- modTargeted1$GeneID %>% unique
reducedExpressionData <- geneExpressionForAnalysis$rosmapDLPFC[,genes]
dim(reducedExpressionData)

pheatmap::pheatmap(cor(reducedExpressionData[,sample(1:ncol(reducedExpressionData),500)]),show_rownames = F, show_colnames = F)

####make pretty heatmaps with module definitions as annotations