#####pull module set

# schema for table:
# 1. ModuleNameFull
# 2. Module
# 3. ModuleMethod
# 4. ModuleBrainRegion
# 5. GeneSetName
# 6. GeneSetCategoryName
# 7. GeneSetAssociationStatistic
# 8. GeneSetEffect
# 9. GeneSetBrainRegion
# 10. GeneSetDirectionAD
# 11. GeneSetADLinked
# 12. GeneSetAdjustedAssociationStatistic
# 13. EvidenceClass

pullReferenceGeneSets <- function(url){
  #pull gene set from Enrichr
  foo <- RCurl::getURL(url)
  
  #split by carriage returns
  fooSplit <- strsplit(foo,"\n")
  
  #unlist
  fooSplit <- unlist(fooSplit)
  
  #split by tab
  fooSplit <- strsplit(fooSplit,"\t")
  
  #name them by the first element
  catNames <- lapply(fooSplit,function(x) x[1])
  names(fooSplit) <- catNames
  
  #drop the first two elements
  fooSplit <- lapply(fooSplit,function(x) x[-c(1,2)])
  
  fxn1 <- function(x){
    library(dplyr)
    strsplit(x,',') %>%
      sapply(function(y) y[1]) %>%
      return
  }
  #remove 1.0 from elements
  fooSplit2 <- lapply(fooSplit,fxn1)
  names(fooSplit2) <- names(fooSplit)
  return(fooSplit2)
}
splitByBrainRegionAdjustPvalue <- function(x){
  #split by brain region and category to adjust the p-values appropriately given the multiple hypothesis testing burden
  brs <- unique(x$ModuleBrainRegion)
  #print(brs)
  cats <- unique(x$GeneSetCategoryName)
  #print(cats)
  combined <- expand.grid(brs,cats)
  fxn1 <- function(y,z,t){
    foo1 <- dplyr::filter(t,ModuleBrainRegion==y & GeneSetCategoryName==z)
    foo1 <- dplyr::mutate(foo1,GeneSetAdjustedAssociationStatistic = p.adjust(GeneSetAssociationStatistic,method='fdr'))
    return(foo1)
  }
  foo2 <- mapply(fxn1,
                 combined[,1],
                 combined[,2],
                 MoreArgs = list(t=x),
                 SIMPLIFY = FALSE)
  foo2 <- do.call(rbind,foo2)
  return(foo2)
}
getAdGenetics <- function(){
  moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  
  
  #magma enrichments
  magmaResults <- synapseClient::synTableQuery("SELECT * FROM syn10380432")@values
  magmaResults <- dplyr::select(magmaResults,SET,BETA,P)
  colnames(magmaResults) <- c('ModuleNameFull',
                              'GeneSetEffect',
                              'GeneSetAssociationStatistic')
  magmaResults$GeneSetName <- rep('MAGMA',nrow(magmaResults))
  magmaResults$GeneSetCategoryName <- rep('genetics',nrow(magmaResults))
  magmaResults$GeneSetADLinked <- rep(TRUE,nrow(magmaResults))
  magmaResults$GeneSetBrainRegion <- rep(NA,nrow(magmaResults))
  magmaResults$GeneSetDirectionAD <- rep(NA,nrow(magmaResults))
  
  magmaResults <- dplyr::select(magmaResults,
                                ModuleNameFull,
                                GeneSetName,
                                GeneSetCategoryName,
                                GeneSetAssociationStatistic,
                                GeneSetEffect,
                                GeneSetBrainRegion,
                                GeneSetDirectionAD,
                                GeneSetADLinked)
  #####merge magma results with the module set given the schema
  moduleSummary <- dplyr::left_join(moduleSet,
                                    magmaResults)
  
  #add adjustd pvalue by brain region and gene set category
  moduleSummary <- splitByBrainRegionAdjustPvalue(moduleSummary)
  
  #####igap, dbgap, and genecards for ad genes enrichments
  source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  dbgap <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=dbGaP")
  
  
  adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
  adList <- c(adList,'HLA-DRB5','HLA-DRB1')
  adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
  adList <- list(igap = adList)
  adList$dbgap <- dbgap$`Alzheimer Disease`
  
  genecardsObj <- synapseClient::synGet('syn10507702')
  genecards <- data.table::fread(synapseClient::getFileLocation(genecardsObj),data.table=F)
  
  adList$genecards <- genecards$`Gene Symbol`
  
  
  
  adTest <- run_amp_ad_enrichment(adList,
                                  'genetics',
                                  manifestId='syn10338156')
  
  adTest <- dplyr::select(adTest,
                          ModuleNameFull,
                          category,
                          geneSet,
                          fisherPval,
                          fisherOR)
  colnames(adTest) <- c('ModuleNameFull',
                        'GeneSetName',
                        'GeneSetCategoryName',
                        'GeneSetAssociationStatistic',
                        'GeneSetEffect')
  
  adTest$GeneSetADLinked <- rep(TRUE,nrow(adTest))
  adTest$GeneSetBrainRegion <- rep(NA,nrow(adTest))
  adTest$GeneSetDirectionAD <- rep(NA,nrow(adTest))
  adTestSummary <- dplyr::left_join(moduleSet,
                                    adTest)
  adTestSummary <- splitByBrainRegionAdjustPvalue(adTestSummary)
  View(adTestSummary)
  
  
  
  moduleSummary <- rbind(moduleSummary,adTestSummary)
  moduleSummary$EvidenceClass <- 'genetics'
  
  return(moduleSummary)
}

getDeg <- function(){
  moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  degResObj <- synapseClient::synGet("syn10496554")
  load(degResObj@filePath)
  source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  degResults <- run_amp_ad_enrichment(amp.ad.de.geneSets,
                                      "degs",
                                      hgnc = TRUE)
  parseDegName <- function(x){
    library(dplyr)
    foo1 <- strsplit(x,'\\.')[[1]]
    br <- foo1[1]
    dir <- foo1[length(foo1)]
    #cate <- foo1[2]
    cate <- paste0(foo1[2:(length(foo1) - 1)],collapse = '_')
    #if(length(grep(paste0('.',br,'_'),cate)) > 0) {
    #  cate <- gsub(paste0('.',br,'_'),'.',cate)
    #}
    
    c('brainRegion' = br,
      'Direction' = dir,
      'reducedCategory' = cate,
      'Category' = x) %>% return
  }
  
  categoryKey <- sapply(unique(degResults$category),
                        parseDegName)
  categoryKey <- t(categoryKey)
  categoryKey <- data.frame(categoryKey,stringsAsFactors = F)
  
  degResults2 <- dplyr::left_join(degResults,categoryKey,by = c('category' = 'Category'))
  degResultsModified <- dplyr::select(degResults2,
                                      ModuleNameFull,
                                      reducedCategory,
                                      geneSet,
                                      fisherPval,
                                      fisherOR,
                                      brainRegion,
                                      Direction)
  
  colnames(degResultsModified) <- c('ModuleNameFull',
                                    'GeneSetName',
                                    'GeneSetCategoryName',
                                    'GeneSetAssociationStatistic',
                                    'GeneSetEffect',
                                    'GeneSetBrainRegion',
                                    'GeneSetDirectionAD')
  
  moduleSummaryDeg <- dplyr::left_join(moduleSet,degResultsModified)
  
  moduleSummaryDeg$GeneSetBrainRegion[moduleSummaryDeg$GeneSetBrainRegion=='CER']<-'CBE'
  
  ###match brain regions for clarity sake
  moduleSummaryDeg <- dplyr::filter(moduleSummaryDeg,ModuleBrainRegion==GeneSetBrainRegion)
  
  ##define which ones are ad related
  #ad_related <- grep('AD',moduleSummaryDeg$GeneSetName)
  moduleSummaryDeg$GeneSetADLinked <- rep(TRUE,nrow(moduleSummaryDeg))
  #moduleSummaryDeg$GeneSetADLinked[ad_related] <- TRUE
  
  moduleSummaryDeg <- splitByBrainRegionAdjustPvalue(moduleSummaryDeg)
  moduleSummaryDeg$EvidenceClass <- rep('deg',nrow(moduleSummaryDeg))
  return(moduleSummaryDeg)
}

getPathways <- function(){}

getEigengene <- function(){}

getModulePreservation <- function(){}

