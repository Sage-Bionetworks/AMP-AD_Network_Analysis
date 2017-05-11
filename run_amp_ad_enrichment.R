#########run enrichment function for manifest
run_amp_ad_enrichment <- function(geneSetList,
                                  geneSetName,
                                  manifestId = "syn9770791"){
  library(dplyr)
  cat('logging into Synapse...\n')
  synapseClient::synapseLogin()
  #grab module definitions
  cat('pulling modules...\n')
  allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ",manifestId))@values
  
  listify <- function(x,y,z){
    ###fxn will listify a long form table
    ###x: unique key
    ###y: values
    ###z: keys
    return(unique(y[which(z==x)]))
  }
  cat('building module gene sets...\n')
  modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                             listify,
                             allMods$external_gene_name,
                             allMods$ModuleNameFull)
  names(modulesLargeList) <- unique(allMods$ModuleNameFull)
  cat('removing genes that are not relevant from reference set...\n')
  #get unique gene keys,drop categories in both cases that are 0 in size
  uniqueModuleList <- modulesLargeList %>%
    unlist %>%
    unique
  
  uniqueGeneSet <- geneSetList %>%
    unlist %>%
    unique
  
  refGeneSet <- uniqueModuleList
  
  cat('running enrichments....\n')
  
  res <- list()
  res$fisher <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapper,
                                            modulesLargeList,
                                            geneSetList,
                                            refGeneSet)
  
  res$pval <- res$fisher[which(1:nrow(res$fisher)%%2==1),]
  rownames(res$pval) <- names(geneSetList)
  res$OR <- res$fisher[which(1:nrow(res$fisher)%%2==0),]
  rownames(res$OR) <- names(geneSetList)
  
  cat('producing tidy data frame....\n')
  
  pval1 <- t(res$pval)
  pval1 <- data.frame(pval1,stringsAsFactors=F)
  pval1 <- dplyr::mutate(pval1,ModuleNameFull = rownames(pval1))
  gatherTest1 <- tidyr::gather(pval1,category,fisherPval,-ModuleNameFull)
  
  or1 <- t(res$OR)
  or1 <- data.frame(or1,stringsAsFactors=F)
  or1 <- dplyr::mutate(or1,ModuleNameFull = rownames(or1))
  gatherTest2 <- tidyr::gather(or1,category,fisherOR,-ModuleNameFull)
  
  gatherTest <- dplyr::left_join(gatherTest1,
                                 gatherTest2)
  gatherTest <- dplyr::mutate(gatherTest,
                              geneSet = geneSetName)
  return(gatherTest)
  
}
  
  
  