#pull down consensus modules
moduleDfList <- rSynapseUtilities::loadDelimIntoList('syn10146717')
moduleDf <- moduleDfList[[1]]
moduleDf <- dplyr::select(moduleDf,-V1)

extractTissueType <- function(brainRegion,
                              fullDf){
  set1 <- grep(brainRegion,fullDf$ModuleNameFull)
  set2 <- grep(brainRegion,fullDf$category)
  set3 <- intersect(set1,set2)
  brDf <- fullDf[set3,]
  return(brDf)
}
dlpfc <- extractTissueType('DLPFC',
                           moduleDf)

dlpfc <- dplyr::filter(dlpfc,p.adjust(dlpfc$fisherPval,method='fdr')<=0.05)
