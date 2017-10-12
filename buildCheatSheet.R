synapseClient::synapseLogin()

source('moduleAnalysis/summaryManifestFunctions.R')

adGeneticsSummary <- getAdGenetics(synId='syn10915669')
adGeneticsSummary <- dplyr::filter(adGeneticsSummary,GeneSetAdjustedAssociationStatistic<=0.05)

degSummary <- getDeg(synId='syn10915669')
degSummary <- dplyr::filter(degSummary,GeneSetAdjustedAssociationStatistic<=0.05)

cellSummary <- getCellTypes(synId='syn10915669')
cellSummary <- dplyr::filter(cellSummary,GeneSetAdjustedAssociationStatistic<=0.05)


adGeneticsCollapsed <- dplyr::select(adGeneticsSummary,
                                     ModuleNameFull,
                                     GeneSetName,
                                     GeneSetADLinked)

degSummary$GeneSetName <- paste0(degSummary$GeneSetName,
                                 degSummary$GeneSetDirectionAD)
degCollapsed <- dplyr::select(degSummary,
                              ModuleNameFull,
                              GeneSetName,
                              GeneSetADLinked)
cellCollapsed <- dplyr::select(cellSummary,
                               ModuleNameFull,
                               GeneSetName,
                               GeneSetADLinked)

fullDf <- rbind(adGeneticsCollapsed,
                degCollapsed,
                cellCollapsed)

masterSheet <- tidyr::spread(fullDf,
                             GeneSetName,
                             GeneSetADLinked)

rownames(masterSheet) <- masterSheet[,1]
masterSheet <- masterSheet[,-1]
masterSheet[!is.na(masterSheet)] <- TRUE
masterSheet[is.na(masterSheet)] <- FALSE

pheatmap::pheatmap(t(data.matrix(masterSheet)))

colnames(masterSheet)
fixDegName <- function(x,br){
  if(x=='AD_CONTROL_AODDOWN'){
    foo <- paste0(br,'.AD.CONTROL.AOD.DOWN')
  } else if (x=='AD_CONTROL_AODUP'){
    foo <- paste0(br,'.AD.CONTROL.AOD.UP')
  } else if (x=='AD_CONTROL_FEMALEDOWN'){
    foo <- paste0(br,'.AD.CONTROL.FEMALE.DOWN')
  } else if (x=='AD_CONTROL_FEMALEUP'){
    foo <- paste0(br,'.AD.CONTROL.FEMALE.UP')
  } else if (x=='AD_CONTROL_MALEDOWN'){
    foo <- paste0(br,'.AD.CONTROL.MALE.DOWN')
  } else if (x=='AD_CONTROL_MALEUP'){
    foo <- paste0(br,'.AD.CONTROL.MALE.UP')
  } else if (x=='AD_CONTROLDOWN'){
    foo <- paste0(br,'.AD.CONTROL.DOWN')
  } else if (x=='AD_CONTROLUP'){
    foo <- paste0(br,'.AD.CONTROL.UP')
  } else if (x=='ApoE1_ApoE0DOWN'){
    foo <- paste0(br,'.ApoE1.ApoE0.DOWN')
  } else if (x=='ApoE1_ApoE0UP'){
    foo <- paste0(br,'.ApoE1.ApoE0.UP')
  } else if (x=='ApoE2_ApoE0DOWN'){
    foo <- paste0(br,'.ApoE2.ApoE0.DOWN')
  } else if (x=='ApoE2_ApoE0UP'){
    foo <- paste0(br,'.ApoE2.ApoE0.UP')
  } else {
    foo <- x
  }
  return(foo)
}
#colnames(masterSheet)[1:12] <- c('')
moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ","syn10915669"))@values
colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
rownames(moduleSet) <- moduleSet$ModuleNameFull
#build manifest for mapply of pull all results
moduleNames <- rownames(masterSheet)
brainRegions <- moduleSet[moduleNames,'ModuleBrainRegion']
annosFull <- lapply(moduleNames,function(mn,masterSheet){
  return(colnames(masterSheet)[which(masterSheet[mn,]!=0)])
}, masterSheet)
annosFull <- mapply(function(x,br){
  return(sapply(x,fixDegName,br))
},annosFull,brainRegions,SIMPLIFY = FALSE)



data4mp<-getDataForModulePainting()


fullManifest <- mapply(pull_all_results,
                       moduleNames,
                       annosFull,
                       MoreArgs = list(geneExpressionForAnalysis = data4mp$geneExpressionForAnalysis,
                                       moduleSet = data4mp$moduleSet,
                                       combined = data4mp$combined),
                       SIMPLIFY = FALSE)
names(fullManifest) <- moduleNames

###tcx blue tcx
save(fullManifest,file='aggregate_module_mainfest.rda')
#annos <- colnames(masterSheet)[which(masterSheet['aggregateTCXblueTCX',]!=0)]


getTopHobs <- function(x){
  return(dplyr::arrange(x$anno,desc(hubs))$GeneID[1:10])
}

topHubs <- sapply(fullManifest,getTopHobs)
topHubs<-names(sort(table(c(topHubs)),decreasing=T))


res1 <- pull_all_results('aggregateTCXblueTCX',
                 c('Zhang.Astrocyte',
                   'Zhang.Endothelial',
                   'Zhang.OPC',
                   'pantherPresenilin',
                   'genecards',
                   'dbgap',
                   'TCX.ApoE1.ApoE0.UP',
                   'TCX.ApoE2.ApoE0.UP',
                   'jensenDisease',
                   'TCX.AD.CONTROL.FEMALE.UP',
                   'TCX.AD.CONTROL.MALE.UP',
                   'TCX.AD.CONTROL.UP'))

str(res1)
#res1$rowAnnoMat$logitDiagnosis[is.na(res1$rowAnnoMat$logitDiagnosis)]<-2

ram <- res1$rowAnnoMat
rownames(ram) <- ram$aSampleId
ram <- dplyr::select(ram,-aSampleId)
rownames(res1$rowAnnoMat) <- res1$rowAnnoMat$aSampleId
res1$rowAnnoMat <- dplyr::select(res1$rowAnnoMat,logitDiagnosis)
png(file='~/Desktop/bluetcx.png',
    height=1600,
    width=2400)

pheatmap::pheatmap(res1$expr,
                   annotation_col = data.frame(data.matrix(res1$colAnno))[,c(1,2,4,6,7)],
                   #annotation_row = data.frame(res1$rowAnnoMat),
                   annotation_row = ram,
                   show_colnames = F, 
                   show_rownames = F)
dev.off()
