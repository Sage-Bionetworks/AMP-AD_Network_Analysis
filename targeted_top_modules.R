####login to Synapse
synapseClient::synapseLogin()

####pull top rankings from table
bar <- synapseClient::synTableQuery("SELECT * FROM syn10516371")@values
View(bar)

####filter out consensus
bar <- dplyr::filter(bar,ModuleMethod!='consensus')

####pull module definitions
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10338156")@values



####pull expression data-set
source('dataPulling/pullExpressionAndPhenoWinsorized.R')

####run consensus clustering per tissue type via interesting modules
tissue <- 'DLPFC'

####pull pairwise relationships between modules
pairwise <- synapseClient::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like \'%DLPFC\' and category like \'%DLPFC\'")@values






keep_mods <- dplyr::filter(bar,ModuleBrainRegion==tissue)
View(keep_mods)
pairwise <- dplyr::filter(pairwise,ModuleNameFull%in%keep_mods$ModuleNameFull & category %in% keep_mods$ModuleNameFull)

pairwise <- utilityFunctions::removeSwappedDupKeyValueDf(pairwise)
pairwise <- dplyr::mutate(pairwise,adj=p.adjust(fisherPval,method='fdr'))
pairwise <- dplyr::filter(pairwise,adj<=0.05)


pairwise <- dplyr::mutate(pairwise, weight = 1/fisherOR)
write.csv(pairwise,file='~/Desktop/dlpfc_targeted.csv',quote=F)


graph1 <- igraph::graph_from_data_frame(pairwise,directed=FALSE)

test1 <- igraph::optimal.community(graph1)

metaGraph <- data.frame(ModuleNameFull = test1$names,
                        metaModule = test1$membership,
                        stringsAsFactors=F)
write.csv(metaGraph,file='~/Desktop/metaGraphDLPFC.csv',quote=F)

modTargeted1 <- dplyr::filter(allMods,ModuleNameFull %in% metaGraph$ModuleNameFull[metaGraph$metaModule==1])

intersectByModule <- function(df){
  uniqueMods <- unique(df$ModuleNameFull)
  geneSet1 <- df$GeneID[df$ModuleNameFull==uniqueMods[1]]
  for ( i in 2:length(uniqueMods)){
    geneSet1 <- intersect(geneSet1,df$GeneID[df$ModuleNameFull==uniqueMods[i]])
  }
  return(geneSet1)
}
keep1a <- intersectByModule(modTargeted1)


View(modTargeted)

library(dplyr)
genes <- modTargeted1$GeneID %>% unique
reducedExpressionData <- geneExpressionForAnalysis$rosmapDLPFC[,genes]
dim(reducedExpressionData)

pheatmap::pheatmap(cor(reducedExpressionData[,sample(1:ncol(reducedExpressionData),500)]),show_rownames = F, show_colnames = F)

####make pretty heatmaps with module definitions as annotations