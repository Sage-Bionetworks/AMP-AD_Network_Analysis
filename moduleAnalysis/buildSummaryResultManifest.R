#####script to build module summary table

synapseClient::synapseLogin()

source('moduleAnalysis/summaryManifestFunctions.R')
#####igap gwas enrichments
adGeneticsSummary <- getAdGenetics()

#####deg enrichments
degSummary <- getDeg()



adGeneticsSummarySig <- dplyr::filter(adGeneticsSummary, 
                                      GeneSetAdjustedAssociationStatistic <= 0.05)
admodcheat <- dplyr::select(adGeneticsSummarySig,
                            ModuleNameFull,
                            GeneSetName,
                            GeneSetADLinked)

admodcheat <- tidyr::spread(admodcheat,
                            ModuleNameFull,
                            GeneSetADLinked)
rownames(admodcheat) <- admodcheat$GeneSetName
admodcheat <- admodcheat[,-1]
admodcheat <- t(admodcheat)
admodcheat[is.na(admodcheat)] <- 0
admodcheat <- data.frame(admodcheat,stringsAsFactors=F)
admodcheat$adGeneticScore <- rowMeans(admodcheat)
admodcheat$ModuleNameFull <- rownames(admodcheat)
admodcheat <- dplyr::arrange(admodcheat,desc(adGeneticScore))

moduleSummarySig <- dplyr::filter(degSummary,
                                  GeneSetAdjustedAssociationStatistic <=0.05)

library(dplyr)
getModuleCheatSheet <- dplyr::select(moduleSummarySig,
                                     ModuleNameFull,
                                     GeneSetName,
                                     GeneSetDirectionAD,
                                     GeneSetBrainRegion,
                                     GeneSetCategoryName,
                                     GeneSetADLinked)
getModuleCheatSheet$genesetdir <- paste0(getModuleCheatSheet$GeneSetName,
                                         getModuleCheatSheet$GeneSetDirectionAD,
                                         getModuleCheatSheet$GeneSetBrainRegion,
                                         getModuleCheatSheet$GeneSetCategoryName)

getModuleCheatSheet <- dplyr::select(getModuleCheatSheet,
                                     ModuleNameFull,
                                     genesetdir,
                                     GeneSetADLinked)

moduleCheatSheet <- tidyr::spread(getModuleCheatSheet,
                                  ModuleNameFull,
                                  GeneSetADLinked)

rownames(moduleCheatSheet) <- moduleCheatSheet$genesetdir
moduleCheatSheet <- moduleCheatSheet[,-1]
moduleCheatSheet <- t(moduleCheatSheet)


moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')

#dropCols <- which(apply(moduleCheatSheet,2,sum,na.rm=T)==0)
#moduleCheatSheet <- moduleCheatSheet[,-dropCols]
moduleCheatSheet[is.na(moduleCheatSheet)] <- 0
moduleCheatSheet <- data.frame(moduleCheatSheet,stringsAsFactors=F)
moduleCheatSheet$degScore <- rowMeans(moduleCheatSheet)
moduleCheatSheet$ModuleNameFull <- rownames(moduleCheatSheet)
moduleCheatSheet <- dplyr::arrange(moduleCheatSheet,desc(degScore))

combinedScores <- dplyr::select(moduleCheatSheet,ModuleNameFull,degScore)
combinedScores$degScore <- as.numeric(scale(combinedScores$degScore,center = FALSE))
combinedScores <- dplyr::full_join(combinedScores,dplyr::select(admodcheat,ModuleNameFull,adGeneticScore))
combinedScores$adGeneticScore <- as.numeric(scale(combinedScores$adGeneticScore,center=FALSE))
combinedScores$aggregate <- combinedScores$degScore + combinedScores$adGeneticScore
combinedScores <- dplyr::arrange(combinedScores,desc(aggregate))
combinedScoresReducted <- combinedScores[1:263,]
combinedScoresReducted <- dplyr::left_join(combinedScoresReducted,moduleSet)

rSynapseUtilities::makeTable(combinedScoresReducted,tableName = "top amp-ad mods september 21 2017",projectId = 'syn5569099')

degTopScores <- dplyr::select(moduleCheatSheet,degScore,ModuleNameFull)
degTopScores <- dplyr::left_join(degTopScores,moduleSet)
rSynapseUtilities::makeTable(degTopScores,tableName = "deg mods september 26 2017",projectId = 'syn5569099')
########compile annotations of modules
####download mods
combinedScoresReducted <- synapseClient::synTableQuery("SELECT * FROM syn10516371")@values


cell_type <- getCellTypes()
cellSummarySig <- dplyr::filter(cell_type, 
                                      GeneSetAdjustedAssociationStatistic <= 0.05)



####get cell type enrichments
#cellTypeEnrichments <- synapseClient::synTableQuery("SELECT * FROM syn10498382")@values

#colnames(cellTypeEnrichments)[c(2,3,7)] <- c("GeneSetCategoryName",
#                                             "GeneSetAssociationStatistic",
#                                             "ModuleBrainRegion")
#cellTypeEnrichments <- splitByBrainRegionAdjustPvalue(cellTypeEnrichments)
#cellTypeEnrichments <- dplyr::filter(cellTypeEnrichments,GeneSetAdjustedAssociationStatistic<=0.05)
combinedScoresReducted2 <- dplyr::left_join(combinedScoresReducted,dplyr::select(cellSummarySig,ModuleNameFull,GeneSetName,GeneSetEffect))



module_ad_score <- apply(moduleCheatSheet,1,sum)
sort(module_ad_score,decreasing=T)[1:30]


#####pathway annotations
pathwaySummary <- getPathways()

#####TO DO
#####eigengene associations
eigengeneSummary <- getEigengene()

#####mod pres
modulePreservationSummary <- getModulePreservation()





#####compile enrichments
enrichments <- synapseClient::synTableQuery("SELECT * FROM syn10492048")@values
enrichments2 <- dplyr::select(enrichments,ModuleNameFull,category,geneSet,fisherPval,fisherOR)
colnames(enrichments2) <- c('ModuleNameFull',
                           'GeneSetName',
                           'GeneSetCategoryName',
                           'GeneSetAssociationStatistic',
                           'GeneSetEffect')

enrichments2$GeneSetBrainRegion <- rep(NA,nrow(enrichments2))
enrichments2$GeneSetDirectionAD <- rep(NA,nrow(enrichments2))
ad_lists <- grep('alzheimer',(enrichments2$GeneSetName))
#ad_lists2 <- grep('load',unique(enrichments2$GeneSetName))
#ad_lists3 <- grep("AD",unique(enrichments2$GeneSetName))
#ad_lists<-grep('',enrichments2$GeneSetName)
enrichments2$GeneSetADLinked <- rep(FALSE,nrow(enrichments2))
enrichments2$GeneSetADLinked[ad_lists] <- TRUE

#bonferroni womp womp

bonferroni_fun <- function(x,ntests=1e8){
  return(min(1,x*ntests))
}

enrichments2$GeneSetAdjustedAssociationStatistic <- sapply(enrichments2$GeneSetAssociationStatistic,bonferroni_fun)

enrichments2 <- dplyr::left_join(moduleSet,enrichments2)
enrichments3 <- dplyr::filter(enrichments2,GeneSetAdjustedAssociationStatistic <= 0.05)

combinedScoresReducted3 <- dplyr::left_join(combinedScoresReducted2,dplyr::select(enrichments3,ModuleNameFull,GeneSetName,GeneSetCategoryName,GeneSetEffect),by='ModuleNameFull')




ad_lists2 <- which(enrichments2$GeneSetADLinked)

moduleSummary <- rbind(moduleSummary,enrichments2)

moduleSummarySig <- dplyr::filter(moduleSummary,GeneSetAssociationStatistic <=0.05)

library(dplyr)
getModuleCheatSheet <- dplyr::select(moduleSummarySig,
                                     ModuleNameFull,
                                     GeneSetName,
                                     GeneSetDirectionAD,
                                     GeneSetBrainRegion,
                                     GeneSetCategoryName,
                                     GeneSetADLinked)
getModuleCheatSheet$genesetdir <- paste0(getModuleCheatSheet$GeneSetName,
                                         getModuleCheatSheet$GeneSetDirectionAD,
                                         getModuleCheatSheet$GeneSetBrainRegion,
                                         getModuleCheatSheet$GeneSetCategoryName)

getModuleCheatSheet <- dplyr::select(getModuleCheatSheet,
                                     ModuleNameFull,
                                     genesetdir,
                                     GeneSetADLinked)

moduleCheatSheet <- tidyr::spread(getModuleCheatSheet,
                                  ModuleNameFull,
                                  GeneSetADLinked)

rownames(moduleCheatSheet) <- moduleCheatSheet$genesetdir
moduleCheatSheet <- moduleCheatSheet[,-1]
moduleCheatSheet <- t(moduleCheatSheet)

dropCols <- which(apply(moduleCheatSheet,2,sum,na.rm=T)==0)
moduleCheatSheet <- moduleCheatSheet[,-dropCols]
module_ad_score <- apply(moduleCheatSheet,1,sum,na.rm=T)
sort(module_ad_score,decreasing=T)[1:30]


enrichmentManifest <- synapseClient::synTableQuery("SELECT * FROM syn10468216")@values

foobar <- readRDS(synapseClient::synGet(enrichmentManifest$id[1])@filePath)
View(enrichmentManifest)



# summaryDegManifest <- dplyr::group_by(degResults2,
#                                       ModuleNameFull,
#                                       Direction,
#                                       reducedCategory) %>% 
#   dplyr::summarise(medianZ = median(Z),
#                    medianOR = median(fisherOR),
#                    medianPval=median(fisherPval))
# 
# 
# summaryDegManifest <- dplyr::mutate(summaryDegManifest,
#                                     adj = p.adjust(medianPval,
#                                                    method='fdr'))
# 
# summaryDegManifest2 <- dplyr::filter(summaryDegManifest,
#                                      adj<=0.05)
# 
# View(summaryDegManifest2)
# 
# g <- ggplot2::ggplot(summaryDegManifest, 
#                      ggplot2::aes(x=Direction,
#                                   y=medianZ,
#                                   fill=reducedCategory))
# g <- g + ggplot2::geom_boxplot(position='dodge')
# #g <- g + ggplot2::scale_y_log10()
# g <- g + ggplot2::theme_grey(base_size = 20) 
# g

#extract 


#dplyr::summarise(numberOfGenes=length(ModuleName)

#categoryKey <- categoryKey[!duplicated(categoryKey),]


#####cell type results
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
cellMarkers <- GeneSets$Cell_Markers
cellTypeResults <- run_amp_ad_enrichment(cellMarkers,
                                        "celltypes",
                                        hgnc=TRUE)
#####combined manifest
fullManifest <- rbind(degResults,
                      cellTypeResults)

magmaReformat <- dplyr::select(magmaResults,SET,P)

colnames(magmaReformat) <- c('ModuleNameFull','magmaPval')
magmaReformat <- dplyr::mutate(magmaReformat,magmaZ=qnorm(magmaPval,lower.tail=F))

summaryDegManifest2 <- dplyr::left_join(summaryDegManifest,magmaReformat)
summaryDegManifest2 <- dplyr::mutate(summaryDegManifest2,combZ=magmaZ/2+medianZ/2)
summaryDegManifest2 <- dplyr::arrange(summaryDegManifest2,desc(medianZ))
#split by each category
fxn1 <- function(x,y){
  foobar <- dplyr::filter(y,reducedCategory==x)
  return(foobar)
}
splitSummaries <- lapply(unique(summaryDegManifest2$reducedCategory),fxn1,summaryDegManifest2)
names(splitSummaries) <- unique(summaryDegManifest2$reducedCategory)
View(splitSummaries[[1]])
splitSummaries2 <- lapply(splitSummaries,function(x){
  x <- dplyr::arrange(x,desc(medianZ))
  xup <- dplyr::filter(x,Direction=='UP')
  xdown <- dplyr::filter(x,Direction=='DOWN')
  return(rbind(xup[1:5,],xdown[1:5,]))
})

View(splitSummaries2[[9]])


##just take top from each up/down
fxn2 <- function(x){
  foobar1 <- dplyr::filter(x,Direction=='DOWN')
  foobar2 <- dplyr::filter(x,Direction=='UP')
  return(c('down_mod'=foobar1$ModuleNameFull[1],
           'up_mod'=foobar2$ModuleNameFull[1]))
}
getMods <- sapply(splitSummaries,
                  fxn2)
topMods <- t(getMods)

#magmaReformat$category <- rep('MAGMA',nrow(magmaReformat))
#magmaReformat$fisherOR <- rep(NA,nrow(magmaReformat))
fullManifest <- dplyr::select(fullManifest,ModuleNameFull,category,fisherPval,fisherOR)
fullManifest <- rbind(fullManifest,magmaReformat)
fullManifest <- dplyr::mutate(fullManifest,Z = qnorm(fisherPval,lower.tail=F))
fullManifestSquare <- dplyr::select(fullManifest,ModuleNameFull,category,Z)
fullManifestSquare <- tidyr::spread(fullManifestSquare,ModuleNameFull,Z)
rownames(fullManifestSquare) <- fullManifestSquare$category
fullManifestSquare <- dplyr::select(fullManifestSquare, -category)
fullManifestSquare <- data.matrix(fullManifestSquare)
fullManifestSquare[!is.finite(fullManifestSquare)] <- NA
fullManifestSquare <- t(fullManifestSquare)
#fullManifestSquare[is.na(fullManifestSquare)] <- 0
foobar <- apply(fullManifestSquare,1,median,na.rm=T)
#####combined score


#####top modules
