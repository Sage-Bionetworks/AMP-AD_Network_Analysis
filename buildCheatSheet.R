synapseClient::synapseLogin()

source('moduleAnalysis/summaryManifestFunctions.R')

adGeneticsSummary <- getAdGenetics(synId='syn10909862')
adGeneticsSummary <- dplyr::filter(adGeneticsSummary,GeneSetAdjustedAssociationStatistic<=0.05)

degSummary <- getDeg(synId='syn10909862')
degSummary <- dplyr::filter(degSummary,GeneSetAdjustedAssociationStatistic<=0.05)

cellSummary <- getCellTypes(synId='syn10909862')
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
annos <- colnames(masterSheet)[which(masterSheet['aggregateDLPFCblueDLPFC',]!=0)]
res1<-pull_all_results('aggregateDLPFCblueDLPFC',
                 c('Zhang.Astrocyte',
                   'Zhang.Endothelial',
                   'pantherPresenilin',
                   'genecards',
                   'DLPFC.AD.CONTROL.FEMALE.UP',
                   'DLPFC.AD.CONTROL.MALE.DOWN'))

str(res1)
png(file='~/Desktop/bluedlpfc.png',height=1600,width=2400)
pheatmap::pheatmap(res1$expr,
                   annotation_col = data.frame(data.matrix(res1$colAnno)),
                   annotation_row = data.frame(res1$rowAnnoMat),
                   show_colnames = F, 
                   show_rownames = F,
                   fontsize = 30)
dev.off()
