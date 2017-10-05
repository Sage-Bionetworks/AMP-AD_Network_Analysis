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




getAdGenetics <- function(synId='syn10338156'){
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
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
  #kegg, wikipathways, biocarta, panther, jensen_diseases, omim_disease, omim_expanded
  kegg <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016")
  wikipathways <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2016")
  jensen_diseases <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Jensen_DISEASES")
  biocarta <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioCarta_2016")
  panther <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")
  omim_disease <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Disease")
  omim_expanded <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Expanded")
  
  adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
  adList <- c(adList,'HLA-DRB5','HLA-DRB1')
  adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
  adList <- list(igap = adList)
  adList$dbgap <- dbgap$`Alzheimer Disease`
  adList$kegg <- kegg$`Alzheimer's disease_Homo sapiens_hsa05010`
  adList$wikipathwaysMouse <- wikipathways$`Alzheimers Disease_Mus musculus_WP2075`
  adList$wikipathwaysHuman <- wikipathways$`Alzheimers Disease_Homo sapiens_WP2059`
  adList$jensenDisease <- jensen_diseases$`Alzheimer's_disease`
  adList$biocarta <- biocarta$`Deregulation of CDK5 in Alzheimers Disease_Homo sapiens_h_p35alzheimersPathway`
  adList$pantherAmyloid <- panther$`Alzheimer disease-amyloid secretase pathway_Homo sapiens_P00003`
  adList$pantherPresenilin <- panther$`Alzheimer disease-presenilin pathway_Homo sapiens_P00004`
  adList$omim <- omim_disease$`alzheimer_disease`
  adList$omimExpanded <- omim_expanded$`alzheimer_disease`
  
  genecardsObj <- synapseClient::synGet('syn10507702')
  genecards <- data.table::fread(synapseClient::getFileLocation(genecardsObj),data.table=F)
  
  adList$genecards <- genecards$`Gene Symbol`
  
  
  
  adTest <- run_amp_ad_enrichment(adList,
                                  'genetics',
                                  manifestId=synId)
  
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

getDeg <- function(synId='syn10338156'){
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  degResObj <- synapseClient::synGet("syn10496554")
  load(degResObj@filePath)
  source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  degResults <- run_amp_ad_enrichment(amp.ad.de.geneSets,
                                      "degs",
                                      hgnc = TRUE,
                                      manifestId = synId)
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

getCellTypes <- function(synId = 'syn10338156'){
  source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  cell_enrich <- run_amp_ad_enrichment(GeneSets$Cell_Markers,
                                  'cell_type',
                                  manifestId=synId) 
  cell_enrich <- dplyr::select(cell_enrich,
                          ModuleNameFull,
                          category,
                          geneSet,
                          fisherPval,
                          fisherOR)
  colnames(cell_enrich) <- c('ModuleNameFull',
                        'GeneSetName',
                        'GeneSetCategoryName',
                        'GeneSetAssociationStatistic',
                        'GeneSetEffect')
  cell_enrich$GeneSetADLinked <- rep(FALSE,nrow(cell_enrich))
  cell_enrich$GeneSetBrainRegion <- rep(NA,nrow(cell_enrich))
  cell_enrich$GeneSetDirectionAD <- rep(NA,nrow(cell_enrich))
  cell_enrichSummary <- dplyr::left_join(moduleSet,
                                    cell_enrich)
  cell_enrichSummary <- splitByBrainRegionAdjustPvalue(cell_enrichSummary)
  #View(cell_enrichSummary)
  return(cell_enrichSummary)
  #moduleSummary <- rbind(moduleSummary,cell_enrichSummary)
  #moduleSummary$EvidenceClass <- 'cell_type'
}

getPathways <- function(){}

getEigengene <- function(){}

getModulePreservation <- function(){}

pull_all_results <- function(moduleName,annos){
  #####pull expression data
  source('dataPulling/pullExpressionAndPhenoWinsorized.R')
  source('pullADgenes.R')
  names(geneExpressionForAnalysis) <- c('TCX',
                                        'CBE',
                                        'DLPFC',
                                        'FP',
                                        'STG',
                                        'PHG',
                                        'IFG')
  #####get module definition
  #synapseClient::synapseLogin()
  foo <- synapseClient::synTableQuery(paste0("select * from syn10884829 where ModuleNameFull =\'",moduleName,"\'"))@values
  moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10884829")@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  modbr <- moduleSet$ModuleBrainRegion[moduleSet$ModuleName==moduleName]
  adgenes <- pullADgenes()
  deggenes <- pullDegGenes()
  celltype <- pullCellTypeGenes()
  combined <- c(adgenes,
                deggenes,
                celltype)
  
  modIn <- function(x,y){
    return(y%in%x)
  }
  
  paintMod <- sapply(combined,modIn,foo$external_gene_name)
  rownames(paintMod) <- foo$GeneID
  paintMod <- data.frame(paintMod,
                         stringsAsFactors = F)
  paintMod <- paintMod[,annos]

  paintMod$external_gene_name <- rownames(paintMod)
  paintMod <- dplyr::left_join(paintMod,foo)
  rownames(paintMod) <- paintMod$GeneID
  res <- list()
  exprMat <- geneExpressionForAnalysis[[modbr]][,paintMod$GeneID]
  rownames(exprMat) <- geneExpressionForAnalysis[[modbr]]$aSampleId
  res$rowAnnoMat <- geneExpressionForAnalysis[[modbr]][,c('aSampleId','logitDiagnosis')]
  annoMat <- data.matrix(paintMod[,1:length(annos)]) %>% data.frame()
  #annoMat <- apply(annoMat,2,as.factor) %>% data.frame()


  res$expr <- exprMat
  #res$anno <- annoMat
  res$anno <- paintMod
  
  #res$colAnno <- res$colAnno[,c(annos,'hubs')]
  
  bicNets <- synapseClient::synTableQuery(paste0("SELECT * FROM syn8681664 where ( (method = \'bic\') and (tissueTypeAbrv = \'",modbr,"\' )  and ( assay = \'RNAseq\'))"))@values
  
  load(synapseClient::synGet(bicNets$id[1])@filePath)
  res$adjacencyMatrix <- as.matrix(bicNetworks$network[paintMod$GeneID,paintMod$GeneID])
  hubs <- data.frame(hubs = rowSums(res$adjacencyMatrix+t(res$adjacencyMatrix)),GeneID=rownames(res$adjacencyMatrix),stringsAsFactors=F)
  res$anno <- dplyr::left_join(res$anno,hubs)
  res$colAnno <- res$anno
  rownames(res$colAnno) <- res$anno$GeneID
  res$colAnno <- res$colAnno[,c(annos,'hubs')]
  return(res)
}


#bicNets <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( method = 'bic' ) AND ( study = 'ROSMAP' ) )")@values



