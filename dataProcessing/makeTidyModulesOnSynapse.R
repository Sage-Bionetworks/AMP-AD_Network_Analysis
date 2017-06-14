synapseClient::synapseLogin()



tidyModules <- function(brainRegion){
  #get consensus modules for a given brain region
  string1 <- paste0("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'moduleIdentification' OR analysisType = 'consensusModuleIdentification') AND ( method = 'kmeans' OR method = 'wina' OR method = 'speakeasy' OR method = 'megena' OR method = 'wgcna' ) AND (tissueOfOrigin = '",brainRegion,"') AND (columnScaled = 'TRUE') )")
  
  #
  foo <- synapseClient::synTableQuery(string1)@values
  bar <- lapply(foo$id,synapseClient::synGet)
  modules <- lapply(bar,function(x){
    library(dplyr)
    synapseClient::getFileLocation(x) %>%
      data.table::fread(data.table=F)})
  
  names(modules) <- foo$method
  modules$speakEasy$V1 <- modules$kmeans$Gene.ID
  colnames(modules$speakEasy) <- c('Gene.ID','Module')
  colnames(modules$kmeans)[1:2] <- c('Gene.ID','Module')
  modules$kmeans <- dplyr::select(modules$kmeans,-moduleLabel)
  
  colnames(modules$megena)[1:2] <- c('Gene.ID','Module')
  modules$wina <- dplyr::select(modules$wina,Geneid,module)
  colnames(modules$wina)[1:2] <- c('Gene.ID','Module')
  tidyDf <- mapply(function(x,y){
    return(cbind(x,rep(y,nrow(x))))
  },
  modules,
  names(modules),
  SIMPLIFY=FALSE)
  tidyDf <- do.call(rbind,tidyDf)
  colnames(tidyDf)[3] <- 'method'
  
  #make custom name
  tidyDf <- dplyr::mutate(tidyDf,ModuleName = paste0(tidyDf$method,tidyDf$Module))
  
  #get hugo ids
  ensgIds <- unique(tidyDf$Gene.ID)
  fullIds <- utilityFunctions::convertEnsemblToHgnc(ensgIds)
  
  tidyDf <- dplyr::left_join(tidyDf,fullIds,by=c('Gene.ID'='ensembl_gene_id'))
  colnames(tidyDf)[1] <- 'GeneID'
  
  return(tidyDf)
}

DLPFC <- tidyModules('dorsolateralPrefrontalCortex')
IFG <- tidyModules('inferiorFrontalGyrus')
STG <- tidyModules('superiorTemporalGyrus')
TCX <- tidyModules('temporalCortex')
CER <- tidyModules('cerebellum')
PHG <- tidyModules('parahippocampalGyrus')
FP <- tidyModules('frontalPole')

rSynapseUtilities::makeTable(IFG,'msbb ifg modules','syn2370594')
rSynapseUtilities::makeTable(STG,'msbb stg modules','syn2370594')
rSynapseUtilities::makeTable(TCX,'mayornaseq tcx modules','syn2370594')
rSynapseUtilities::makeTable(CER,'mayornaseq cer modules','syn2370594')
rSynapseUtilities::makeTable(PHG,'msbb phg modules','syn2370594')
rSynapseUtilities::makeTable(FP,'msbb fp modules','syn2370594')
