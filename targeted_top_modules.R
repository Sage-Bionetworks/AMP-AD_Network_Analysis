####login to Synapse
synapseClient::synapseLogin()

####pull top rankings from table
buildTargetedModules <- function(tissueType,synId="syn10909307"){
  library(dplyr)
  bar <- synapseClient::synTableQuery(paste0("SELECT * FROM ",synId))@values %>%
    dplyr::filter(ModuleMethod!='consensus' & ModuleBrainRegion==tissueType)
  moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
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
  res$moduleGraph$weight <- res$moduleGraph$weight+1e-15
  graph1 <- igraph::graph_from_data_frame(res$moduleGraph,directed=FALSE)
  set.seed(2)
  test1 <- igraph::edge.betweenness.community(graph1)
  metaGraph <- data.frame(ModuleNameFull = test1$names,
                          metaModule = test1$membership,
                          stringsAsFactors=F)
  metaGraph2<-dplyr::left_join(metaGraph,moduleSet)
  
  
  bootstrapFrequencies <- function(i,df,allMods,tissueType){
    set.seed(47)
    print(i)
    bar1 <- dplyr::filter(df,metaModule==i)

    bar2 <- lapply(unique(df$ModuleMethod),function(x,y){
      library(dplyr)
      dplyr::filter(y,
                    ModuleMethod==x) %>%
      return()
    },df)
    names(bar2) <- unique(df$ModuleMethod)
    bar3 <- table(bar1$ModuleMethod)
    bar3 <- bar3[names(bar2)]
    if(sum(is.na(bar3))>0){
      bar2 <- bar2[!is.na(bar3)]
      bar3 <- bar3[!is.na(bar3)]
    }
    fxn1 <- function(modCount,
                     df){
      foo1 <- sapply(1:100,
                     function(x,dta,modcount){
        ind <- sample(1:nrow(dta),modcount)
        return(dta$ModuleNameFull[ind])
      },df,modCount)
      foo1 <- as.matrix(foo1)
      if(nrow(foo1)==100){
        foo1 <- t(foo1)
      }
      return(foo1)
      #sapply(1:1000,sample(1:nrow(df),modCount))
    }
    bar4 <- mapply(fxn1,
                   bar3,
                   bar2,
                   SIMPLIFY=F)
    getFrequencies <- function(modNameMat,allMods,refGene){
      library(dplyr)
      apply(modNameMat,2,function(x,allMods,refGene){
        genes<-dplyr::filter(allMods,ModuleNameFull %in% x)
        genes <- table(genes$GeneID)
        refGene[names(genes)]<-genes
        return(refGene)
      },allMods,refGene) %>% return
    }
    
    refGene <- unique(dplyr::filter(allMods,brainRegion==tissueType)$GeneID)
    refGene2 <- rep(0,length(refGene))
    names(refGene2) <- refGene
    
    modGenes <- dplyr::filter(allMods,ModuleNameFull %in% bar1$ModuleNameFull)$GeneID
    modGenes2 <- table(modGenes)
    refGene3 <- refGene2
    refGene3[names(modGenes2)] <- modGenes2
    
    
    frequencyMatrices <- lapply(bar4,getFrequencies,allMods,refGene2)
    names(frequencyMatrices) <- names(bar4)
    masterFrequencyMatrix <- Reduce("+",frequencyMatrices)
    fivePercent <- apply(masterFrequencyMatrix,1,quantile,.95)
    gen <- list()
    gen$ensg <- names(which(refGene3>fivePercent))
    return(gen)
  }

  res$mods<-lapply(unique(metaGraph2$metaModule),
         bootstrapFrequencies,
         metaGraph2,
         allMods,
         tissueType)
  
  res$moduleGraphCommunities <- metaGraph2
  
  modLen<-sapply(res$mods,function(x) length(x$ensg))
  res$mods <- res$mods[modLen>0]
  
  # getMajority <- function(df){
  #   masterTableHGNC <- table(df$external_gene_name)
  #   masterTableENSG <- table(df$GeneID)
  #   #masterTable <- masterTable/len1
  #   gen <- list()
  #   masterTableHGNC <- sort(masterTableHGNC,decreasing=T)
  #   masterTableENSG <- sort(masterTableENSG,decreasing=T)
  #   if(length(masterTableHGNC)< 1000 | length(masterTableENSG) < 1000){
  #     gen$ensg <- names(masterTableENSG)
  #     gen$hgnc <- names(masterTableHGNC)
  #   }else{
  #     gen$ensg <- names(masterTableENSG)[1:1000]
  #     gen$hgnc <- names(masterTableHGNC)[1:1000]
  #   }
  #   #masterTableHGNC <- which(masterTableHGNC > 1)
  #   #masterTableENSG <- which(masterTableENSG > 1)
  # 
  #   return(gen)
  # }
  #38 41 17 53  9 16 27  4 
  mods <- unique(res$moduleGraphCommunities$metaModule)[modLen>0]
  # dfList <- lapply(mods,function(x,df,allMods){
  #   library(dplyr)
  #   dplyr::filter(allMods,ModuleNameFull %in% df$ModuleNameFull[df$metaModule==x]) %>%
  #     return},
  #   res$moduleGraphCommunities,
  #   allMods)
  # 
  # res$mods <- lapply(dfList,getMajority)
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

#rSynapseUtilities::makeTable(AggregateModuleManifest,'collapsed ad meta modules september 21 2017',projectId='syn5569099')

rSynapseUtilities::makeTable(AggregateModuleManifest,'collapsed ad meta modules october 10 2017',projectId='syn2370594')


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