pairwiseComp <- data.table::fread('~/Desktop/pairwisemods.csv',
                                  data.table=F)
pairwiseComp <- dplyr::select(pairwiseComp,-V1)
gc()

#remove self tests
whichSelf <- which(pairwiseComp$ModuleNameFull!=pairwiseComp$category)
pairwiseComp <- pairwiseComp[whichSelf,]

pairwiseComp <- dplyr::mutate(pairwiseComp,
                              adj.pval = p.adjust(fisherPval,method = 'fdr'))
pairwiseCompFilt <- dplyr::filter(pairwiseComp,adj.pval<=0.05)

pairwiseCompFilt <- dplyr::select(pairwiseCompFilt,ModuleNameFull,fisherOR,category)
write.csv(pairwiseCompFilt,
          file='~/Desktop/pairwiseModuleComparisons.csv',
          quote=F)

synapseClient::synapseLogin()
nodeTable <-synapseClient::synTableQuery("SELECT * FROM syn9770842")@values
nodeTable <- dplyr::select(nodeTable,ModuleNameFull,method,brainRegion,adGeneticEnrich,microgEnrich,neurEnrich,astroEnrich,endoEnrich)
nodeTable <- nodeTable[which(!duplicated(nodeTable)),]
write.csv(nodeTable,file='~/Desktop/nodeDescriptions.csv')
library(dplyr)
mostLikelyCell <- dplyr::select(nodeTable,microgEnrich,neurEnrich,astroEnrich,endoEnrich) %>%
  apply(1,function(x) if(min(x)<1e-6){return(which.min(x))}else{return(NA)})
nodeTable2 <- dplyr::mutate(nodeTable,cell=c('microglia','neuron','astrocyte','endothelium')[mostLikelyCell])
write.csv(dplyr::select(nodeTable2,ModuleNameFull,cell),file='~/Desktop/moduleCells.csv',quote=F,row.names=F)
