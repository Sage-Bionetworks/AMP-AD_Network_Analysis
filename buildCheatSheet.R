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


###tcx blue tcx
annos <- colnames(masterSheet)[which(masterSheet['aggregateTCXblueTCX',]!=0)]
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
