#####script to build module summary table

synapseClient::synapseLogin()

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


moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')

#####MAGMA results
magmaResults <- synapseClient::synTableQuery("SELECT * FROM syn10380432")@values
View(magmaResults)
magmaResults <- dplyr::select(magmaResults,SET,BETA,P)
colnames(magmaResults) <- c('ModuleNameFull',
                            'GeneSetEffect',
                            'GeneSetAssociationStatistic')
magmaResults$GeneSetName <- rep('MAGMA',
                                nrow(magmaResults))
magmaResults$GeneSetCategoryName <- rep('genetics',
                                        nrow(magmaResults))
magmaResults$GeneSetADLinked <- rep(TRUE,
                                    nrow(magmaResults))
magmaResults$GeneSetBrainRegion <- rep(NA,
                                       nrow(magmaResults))

magmaResults$GeneSetDirectionAD <- rep(NA,
                                       nrow(magmaResults))

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

View(moduleSummary)
#add adjustd pvalue by brain region and gene set category
splitByBrainRegionAdjustPvalue <- function(x){
  #split by brain region and category to adjust the p-values appropriately given the multiple hypothesis testing burden
  brs <- unique(x$ModuleBrainRegion)
  cats <- unique(x$GeneSetCategoryName)
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

moduleSummary <- splitByBrainRegionAdjustPvalue(moduleSummary)



#####grab degs
#MAYO: syn8468023
#MSSM: syn10157628
#ROSMAP: syn8456721

# mayoResObj <- synapseClient::synGet("syn8468023")
# mayoRes <- data.table::fread(mayoResObj@filePath,
#                                 data.table=FALSE)
# mayoRes <- dplyr::filter(mayoRes,
#                          adj.P.Val <= 0.05)
# mayoRes <- dplyr::mutate(mayoRes,Comparison2 = paste0(Comparison,'_',Direction))
# View(mayoRes)


degResObj <- synapseClient::synGet("syn10163525")
degRes <- data.table::fread(degResObj@filePath,
                            data.table = FALSE)
degRes <- dplyr::filter(degRes,
                        adj.P.Val <= 0.05)
degRes <- dplyr::mutate(degRes,Comparison2 = paste0(Comparison,'_',Direction))

source('enrichmentAnalysis/run_amp_ad_enrichment.R')

listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z == x)]))
}

degList <- lapply(unique(degRes$Comparison2),
                           listify,
                           degRes$Gene.ID,
                           degRes$Comparison2)
names(degList) <- unique(degRes$Comparison2)
reformatNames <- function(x){
  foo1 <- strsplit(x,'\\.')
  foo2 <- sapply(foo1,function(y) y[1])
  return(foo2)
}
degList <- lapply(degList,reformatNames)

degResults <- run_amp_ad_enrichment(degList,
                                    "degs",
                                    hgnc = FALSE)

###reorganize deg results
#split off brain region
parseDegName <- function(x){
  library(dplyr)
  foo1 <- strsplit(x,'_')[[1]]
  br <- foo1[1]
  dir <- foo1[length(foo1)]
  cate <- paste0(foo1[2:(length(foo1) - 1)],collapse = '_')
  if (length(grep(paste0('.',br,'_'),cate)) > 0) {
    cate <- gsub(paste0('.',br,'_'),'.',cate)
  }
  
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

degResults2 <- dplyr::mutate(degResults2,Z = qnorm(fisherPval,lower.tail = F))
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

###match brain regions for clarity sake
moduleSummaryDeg <- dplyr::filter(moduleSummaryDeg,ModuleBrainRegion==GeneSetBrainRegion)

##define which ones are ad related
ad_related <- grep('AD',moduleSummaryDeg$GeneSetName)
moduleSummaryDeg$GeneSetADLinked <- rep(FALSE,nrow(moduleSummaryDeg))
moduleSummaryDeg$GeneSetADLinked[ad_related] <- TRUE

moduleSummaryDeg <- splitByBrainRegionAdjustPvalue(moduleSummaryDeg)

moduleSummary <- rbind(moduleSummary,
                       moduleSummaryDeg)

View(dplyr::filter(moduleSummary,GeneSetAdjustedAssociationStatistic<0.05))


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
