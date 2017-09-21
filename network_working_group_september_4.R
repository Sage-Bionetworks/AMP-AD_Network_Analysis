#login to Synapse
synapseClient::synapseLogin()

#pull enrichments
enrichments <- synapseClient::synTableQuery("SELECT * FROM syn10492048")@values
enrichments2 <- dplyr::select(enrichments,ModuleNameFull,category,geneSet,fisherPval,fisherOR)
colnames(enrichments2) <- c('ModuleNameFull',
                            'GeneSetName',
                            'GeneSetCategoryName',
                            'GeneSetAssociationStatistic',
                            'GeneSetEffect')

enrichments2$GeneSetBrainRegion <- rep(NA,nrow(enrichments2))
enrichments2$GeneSetDirectionAD <- rep(NA,nrow(enrichments2))
ad_lists <- grep('alzheimer',(enrichments2$GeneSetName),ignore.case=TRUE)
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
enrichments3 <- dplyr::filter(enrichments2,GeneSetAdjustedAssociationStatistic <= 0.05)
View(enrichments3)

boo3 <- dplyr::group_by(enrichments3,ModuleNameFull,GeneSetCategoryName)
boo4 <- dplyr::summarise(boo3,mean_or=mean(GeneSetEffect),mean_pval=mean(-log10(GeneSetAssociationStatistic)),nsigcat=length(GeneSetCategoryName))

### get module sizes
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10338156")@values
library(dplyr)
fooSummarize <- dplyr::group_by(allMods,brainRegion,ModuleName,method)%>%
  dplyr::summarise(numberOfGenes=length(ModuleName))

fooSummarize <- dplyr::mutate(fooSummarize,ModuleNameFull = paste0(ModuleName,brainRegion))
module_score_manifest <- dplyr::left_join(boo4,fooSummarize)

module_score_manifest$mean_or[!is.finite(module_score_manifest$mean_or)] <- NA
module_score_manifest<-na.omit(module_score_manifest)
module_score_manifest$adjusted_log10_mean_or <- MASS::rlm(log10(module_score_manifest$mean_or) ~ log10(module_score_manifest$numberOfGenes))$resid
module_score_manifest$adjusted_log10_nsigset <- MASS::rlm(log10(module_score_manifest$nsigcat) ~ log10(module_score_manifest$numberOfGenes))$resid
module_score_manifest$adjusted_mean_pval <- MASS::rlm(module_score_manifest$mean_pval ~ log10(module_score_manifest$numberOfGenes))$resid

module_score_manifest$adjusted_mean_or <- 10^(MASS::rlm(log10(module_score_manifest$mean_or) ~ log10(module_score_manifest$numberOfGenes))$resid)
module_score_manifest$adjusted_nsigset <- 10^(MASS::rlm(log10(module_score_manifest$nsigcat) ~ log10(module_score_manifest$numberOfGenes))$resid)

module_score_manifest <- dplyr::mutate(module_score_manifest,log10_mean_or = log10(mean_or))
module_score_manifest <- dplyr::mutate(module_score_manifest,log10_nsigcat = log10(nsigcat))
#####box plots with dots
#par(mfcol=c(2,1))
makeBoxPlot1 <- function(msm,
                        yaxis,
                        yaxisLabel){
  #str <- paste0("y=log10(",yaxis,")")
  g <- ggplot2::ggplot(msm, 
                       ggplot2::aes_string(x="method",
                                    y=yaxis,
                                    color="method",
                                    fill="method"))
  g <- g + ggplot2::geom_boxplot(alpha = 0.2)
  g <- g + ggplot2::geom_jitter(ggplot2::aes_string(x="method",
                                             y=yaxis,
                                             color="method",
                                             fill="method"))
  #g <- g + ggplot2::scale_y_log10()
  g <- g + ggplot2::theme_grey(base_size = 20) 
  g <- g + ggplot2::coord_flip()
  g <- g + ggplot2::ylab(yaxisLabel)

  return(g)
}
makeBoxPlot1(module_score_manifest,"log10_mean_or", expression(log[10]~average~odds~ratio))
#makeBoxPlot1(module_score_manifest,"adjusted_log10_mean_or")

makeBoxPlot1(module_score_manifest,"log10_nsigcat", expression(log[10]~number~of~significant~pathways))
#makeBoxPlot(module_score_manifest,"adjusted_log10_nsigset")

makeBoxPlot1(module_score_manifest,"mean_pval", expression(average~log[10]~p-value))
#makeBoxPlot(module_score_manifest,"adjusted_mean_pval")


###compute mann whitney u test for all pairs
#extract scores
mean_or<-lapply(unique(module_score_manifest$method),
                utilityFunctions::listify,
                module_score_manifest$log10_mean_or,
                module_score_manifest$method)
names(mean_or) <- unique(module_score_manifest$method)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)

ncat<-lapply(unique(module_score_manifest$method),
                utilityFunctions::listify,
                module_score_manifest$log10_nsigcat,
                module_score_manifest$method)
names(ncat) <- unique(module_score_manifest$method)
testDf<-expand.grid(names(ncat),names(ncat))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=ncat),
  SIMPLIFY=TRUE)

#make Scatter plot
g <- ggplot2::ggplot(module_score_manifest,
                     ggplot2::aes(log10(numberOfGenes),
                                  log10_mean_or))
g <- g + ggplot2::geom_point()
g <- g + ggplot2::geom_smooth()
g <- g + ggplot2::xlab(expression(log[10]~module~size))
g <- g + ggplot2::ylab(expression(log[10]~average_odds_ratio))
g

#after adjustment
g <- ggplot2::ggplot(module_score_manifest,
                     ggplot2::aes(log10(numberOfGenes),
                                  adjusted_log10_mean_or))
g <- g + ggplot2::geom_point()
g <- g + ggplot2::geom_smooth()
g <- g + ggplot2::xlab(expression(log[10]~module~size))
g <- g + ggplot2::ylab(expression(adj~log[10]~average_odds_ratio))
g


g <- ggplot2::ggplot(module_score_manifest,
                     ggplot2::aes(log10(numberOfGenes),
                                  log10_nsigcat))
g <- g + ggplot2::geom_point()
g <- g + ggplot2::geom_smooth()
g <- g + ggplot2::xlab(expression(log[10]~module~size))
g <- g + ggplot2::ylab(expression(log[10]~number~of~significant~pathways))
g

#plot(log10(module_score_manifest$numberOfGenes),module_score_manifest$log10_mean_or)
makeBoxPlot1(module_score_manifest,"adjusted_log10_mean_or", expression(adj~log[10]~average~odds~ratio))

mean_or<-lapply(unique(module_score_manifest$method),
                utilityFunctions::listify,
                module_score_manifest$adjusted_log10_mean_or,
                module_score_manifest$method)
names(mean_or) <- unique(module_score_manifest$method)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)

makeBoxPlot1(module_score_manifest,"adjusted_log10_nsigset", expression(adj~log[10]~number~of~significant~pathways))
ncat<-lapply(unique(module_score_manifest$method),
             utilityFunctions::listify,
             module_score_manifest$adjusted_log10_nsigset,
             module_score_manifest$method)
names(ncat) <- unique(module_score_manifest$method)
testDf<-expand.grid(names(ncat),names(ncat))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=ncat),
  SIMPLIFY=TRUE)



#####fix msm
msm_new <- dplyr::select(module_score_manifest,
                         ModuleNameFull,
                         method,
                         GeneSetCategoryName,
                         log10_mean_or,
                         mean_pval,
                         log10_nsigcat,
                         adjusted_log10_mean_or,
                         adjusted_mean_pval,
                         adjusted_log10_nsigset)

msm_new2 <- tidyr::gather(msm_new,
                          key_or,
                          value_or,
                          log10_mean_or,
                          adjusted_log10_mean_or)

msm_new2 <- tidyr::gather(msm_new2,
                          key_pval,
                          value_pval,
                          mean_pval,
                          adjusted_mean_pval)

msm_new2 <- tidyr::gather(msm_new2,
                          key_cat,
                          value_cat,
                          log10_nsigcat,
                          adjusted_log10_nsigset)
makeBoxPlot <- function(msm,
                        yaxis,
                        key,
                        ylab){
  g <- ggplot2::ggplot(msm, 
                       ggplot2::aes_string(x="method",
                                           y=yaxis,
                                           color=key,
                                           fill=key))
  g <- g + ggplot2::geom_boxplot(alpha = 0.2,
                                 position = 'dodge')
  # g <- g + ggplot2::geom_jitter(ggplot2::aes_string(x="method",
  #                                                    y="value_cat",
  #                                                    color="key_cat",
  #                                                    fill="key_cat"))
  
  #g <- g + ggplot2::scale_y_log10()
  g <- g + ggplot2::theme_grey(base_size = 20) 
  g <- g + ggplot2::coord_flip()
  g <- g + ggplot2::ylab(ylab)
  g
}
makeBoxPlot(msm_new2,"value_or","key_or",expression(log[10]~average~odds~ratio))



##########do same for DEG enrichment
source('moduleAnalysis/summaryManifestFunctions.R')
#####igap gwas enrichments
adGeneticsSummary <- getAdGenetics()

#####deg enrichments
degSummary <- getDeg()

#####do same adjustments for things
###adjustments
moduleSummarySig <- dplyr::filter(degSummary,
                                  GeneSetAdjustedAssociationStatistic <=0.05)

moduleSummarySig <- dplyr::left_join(moduleSummarySig,fooSummarize)
moduleSummarySig <- dplyr::mutate(moduleSummarySig,log10or=log10(GeneSetEffect))
makeBoxPlot1(moduleSummarySig,"log10or",expression(DEG~log[10]~OR))

mean_or<-lapply(unique(moduleSummarySig$ModuleMethod),
                utilityFunctions::listify,
                moduleSummarySig$log10or,
                moduleSummarySig$method)
names(mean_or) <- unique(moduleSummarySig$ModuleMethod)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)

moduleSummarySig$log10or[!is.finite(moduleSummarySig$log10or)] <- NA
moduleSummarySig<-na.omit(moduleSummarySig)
moduleSummarySig$adjusted_log10_mean_or <- MASS::rlm(moduleSummarySig$log10or ~ log10(moduleSummarySig$numberOfGenes))$resid
makeBoxPlot1(moduleSummarySig,"adjusted_log10_mean_or",expression(adj~DEG~log[10]~OR))
mean_or<-lapply(unique(moduleSummarySig$ModuleMethod),
                utilityFunctions::listify,
                moduleSummarySig$adjusted_log10_mean_or,
                moduleSummarySig$method)
names(mean_or) <- unique(moduleSummarySig$ModuleMethod)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)

testDf <- dplyr::arrange(testDf,pval)
testDf <- testDf[1:15,]
testDf <- dplyr::mutate(testDf,adj = p.adjust(pval,method='fdr'))
View(testDf)


#############ad genetics
moduleSummarySig <- dplyr::filter(adGeneticsSummary,
                                  GeneSetAdjustedAssociationStatistic <=0.05)

moduleSummarySig <- dplyr::left_join(moduleSummarySig,fooSummarize)
moduleSummarySig <- dplyr::mutate(moduleSummarySig,log10or=log10(GeneSetEffect))
makeBoxPlot1(moduleSummarySig,"log10or",expression(AD~Genes~log[10]~OR))

mean_or<-lapply(unique(moduleSummarySig$ModuleMethod),
                utilityFunctions::listify,
                moduleSummarySig$log10or,
                moduleSummarySig$method)
names(mean_or) <- unique(moduleSummarySig$ModuleMethod)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)
testDf <- dplyr::arrange(testDf,pval)
testDf <- testDf[1:15,]
testDf <- dplyr::mutate(testDf,adj = p.adjust(pval,method='fdr'))
View(testDf)



moduleSummarySig$log10or[!is.finite(moduleSummarySig$log10or)] <- NA
#moduleSummarySig<-na.omit(moduleSummarySig,"log10or")
moduleSummarySig <- moduleSummarySig[which(!is.na(moduleSummarySig$log10or)),]
moduleSummarySig$adjusted_log10_mean_or <- MASS::rlm(moduleSummarySig$log10or ~ log10(moduleSummarySig$numberOfGenes))$resid
makeBoxPlot1(moduleSummarySig,"adjusted_log10_mean_or",expression(adj~DEG~log[10]~OR))
mean_or<-lapply(unique(moduleSummarySig$ModuleMethod),
                utilityFunctions::listify,
                moduleSummarySig$adjusted_log10_mean_or,
                moduleSummarySig$method)
names(mean_or) <- unique(moduleSummarySig$ModuleMethod)
testDf<-expand.grid(names(mean_or),names(mean_or))
testDf$pval<-mapply(function(x,y,distList){
  return(wilcox.test(distList[[x]],
                     distList[[y]],
                     alternative='greater')$p.value)},
  testDf$Var1,
  testDf$Var2,
  MoreArgs = list(distList=mean_or),
  SIMPLIFY=TRUE)

testDf <- dplyr::arrange(testDf,pval)
testDf <- testDf[1:15,]
testDf <- dplyr::mutate(testDf,adj = p.adjust(pval,method='fdr'))
View(testDf)

#rebuild combined scores with log10 OR instead of just a simple indicator variable

adGeneticsSummarySig <- dplyr::filter(adGeneticsSummary, 
                                      GeneSetAdjustedAssociationStatistic <= 0.05)
admodcheat <- dplyr::select(adGeneticsSummarySig,
                            ModuleNameFull,
                            GeneSetName,
                            GeneSetEffect)

winf <- which(admodcheat$GeneSetEffect==Inf)
maxNonInf <- max(admodcheat$GeneSetEffect[-winf])
admodcheat$GeneSetEffect[winf] <- maxNonInf
admodcheat$GeneSetEffect <- log10(admodcheat$GeneSetEffect)
#which(admodcheat$GeneSetEffect==Inf)

admodcheat <- tidyr::spread(admodcheat,
                            ModuleNameFull,
                            GeneSetEffect)
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
                                     GeneSetEffect)
getModuleCheatSheet$genesetdir <- paste0(getModuleCheatSheet$GeneSetName,
                                         getModuleCheatSheet$GeneSetDirectionAD,
                                         getModuleCheatSheet$GeneSetBrainRegion,
                                         getModuleCheatSheet$GeneSetCategoryName)


getModuleCheatSheet <- dplyr::select(getModuleCheatSheet,
                                     ModuleNameFull,
                                     genesetdir,
                                     GeneSetEffect)

winf <- which(getModuleCheatSheet$GeneSetEffect==Inf)
maxNonInf <- max(getModuleCheatSheet$GeneSetEffect[-winf])
getModuleCheatSheet$GeneSetEffect[winf] <- maxNonInf
getModuleCheatSheet$GeneSetEffect <- log10(getModuleCheatSheet$GeneSetEffect)

moduleCheatSheet <- tidyr::spread(getModuleCheatSheet,
                                  ModuleNameFull,
                                  GeneSetEffect)

rownames(moduleCheatSheet) <- moduleCheatSheet$genesetdir
moduleCheatSheet <- moduleCheatSheet[,-1]
moduleCheatSheet <- t(moduleCheatSheet)




#dropCols <- which(apply(moduleCheatSheet,2,sum,na.rm=T)==0)
#moduleCheatSheet <- moduleCheatSheet[,-dropCols]
moduleCheatSheet[is.na(moduleCheatSheet)] <- 0
moduleCheatSheet <- data.frame(moduleCheatSheet,stringsAsFactors=F)
moduleCheatSheet$degScore <- rowMeans(moduleCheatSheet)
moduleCheatSheet$ModuleNameFull <- rownames(moduleCheatSheet)
moduleCheatSheet <- dplyr::arrange(moduleCheatSheet,desc(degScore))

combinedScores <- dplyr::select(moduleCheatSheet,ModuleNameFull,degScore)
combinedScores$degScore <- as.numeric(scale(combinedScores$degScore))
combinedScores <- dplyr::full_join(combinedScores,dplyr::select(admodcheat,ModuleNameFull,adGeneticScore))
combinedScores$adGeneticScore <- as.numeric(scale(combinedScores$adGeneticScore))
combinedScores$aggregate <- combinedScores$degScore + combinedScores$adGeneticScore
combinedScores <- dplyr::arrange(combinedScores,desc(aggregate))
combinedScoresReducted <- combinedScores[1:233,]
moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")@values
colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')


combinedScoresReducted <- dplyr::left_join(combinedScoresReducted,moduleSet)






#combinedScoresReducted <- synapseClient::synTableQuery("SELECT * FROM syn10516371")@values
combinedScoresReducted <- dplyr::left_join(combinedScoresReducted,fooSummarize)
#View(combinedScoresReducted)

makeBoxPlot1(combinedScoresReducted,"degScore","Scaled Number of Sig DEG enrichments")
makeBoxPlot1(combinedScoresReducted,"adGeneticScore","Scaled Number of Sig AD enrichments")
makeBoxPlot1(combinedScoresReducted,"aggregate","Aggregate AD + DEG enrichments")

cell_type <- getCellTypes()
cellSummarySig <- dplyr::filter(cell_type, 
                                GeneSetAdjustedAssociationStatistic <= 0.05)
combinedScoresReducted2 <- dplyr::left_join(combinedScoresReducted,dplyr::select(cellSummarySig,ModuleNameFull,GeneSetName,GeneSetEffect))


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


#####FOR DLPFC add 
#1) pairwise module relationships
#2) enrichments
#3) cell types

#output to file that can be loaded into cytoscapely
csr_dlpfc <- dplyr::filter(combinedScoresReducted3,brainRegion=='STG')

dlpfc_pairwise <- synapseClient::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%STG' and category like '%STG'")@values
dlpfc_pairwise <- utilityFunctions::removeSwappedDupKeyValueDf(dlpfc_pairwise)
dlpfc_pairwise <- dplyr::mutate(dlpfc_pairwise,adj=p.adjust(fisherPval,method='fdr'))
dlpfc_pairwise <- dplyr::filter(dlpfc_pairwise,adj<=0.05)

uniqueGenes <- unique(csr_dlpfc$ModuleNameFull)
dlpfc_pairwise2 <- dplyr::filter(dlpfc_pairwise,(from%in%uniqueGenes) & (to%in%uniqueGenes) )
dlpfc_pairwise2 <- dplyr::select(dlpfc_pairwise2,from,to,fisherOR)

csr_dlpfc_pathways <- dplyr::select(csr_dlpfc,ModuleNameFull,GeneSetName.y,GeneSetEffect.y)
csr_dlpfc_pathways <- csr_dlpfc_pathways[!duplicated(csr_dlpfc_pathways),]

csr_dlpfc_pathway_anno <- data.frame(unique(csr_dlpfc_pathways$GeneSetName.y),rep('pathway',length(unique(csr_dlpfc_pathways$GeneSetName.y))),stringsAsFactors=F)


csr_dlpfc_cells <- dplyr::select(csr_dlpfc,ModuleNameFull,GeneSetName.x,GeneSetEffect.x)
csr_dlpfc_cells <- csr_dlpfc_cells[!duplicated(csr_dlpfc_cells),]


csr_dlpfc_cells_anno <- data.frame(unique(csr_dlpfc_cells$GeneSetName.x),rep('celltype',length(unique(csr_dlpfc_cells$GeneSetName.x))),stringsAsFactors=F)

csr_dlpfc_deg <- dplyr::filter(moduleSummarySig,ModuleBrainRegion=='STG')
csr_dlpfc_deg$GeneSetName2 <- paste0(csr_dlpfc_deg$GeneSetName,csr_dlpfc_deg$GeneSetDirectionAD)
csr_dlpfc_deg <- dplyr::select(csr_dlpfc_deg,ModuleNameFull,GeneSetName2,GeneSetEffect)

csr_dlpfc_deg_anno <- data.frame(unique(csr_dlpfc_deg$GeneSetName2),rep('DEG',length(unique(csr_dlpfc_deg$GeneSetName2))),stringsAsFactors=F)

csr_dlpfc_ad <- dplyr::filter(adGeneticsSummarySig,ModuleBrainRegion=='STG')
csr_dlpfc_ad <- dplyr::select(csr_dlpfc_ad,ModuleNameFull,GeneSetName,GeneSetEffect)
csr_dlpfc_ad_anno <- data.frame(unique(csr_dlpfc_ad$GeneSetName),rep('AD',length(unique(csr_dlpfc_ad$GeneSetName))),stringsAsFactors=F)


moduleNameAnno <- dplyr::select(moduleSet,ModuleNameFull,ModuleMethod)


masterAnno <- rbind(as.matrix(moduleNameAnno),
                    as.matrix(csr_dlpfc_ad_anno),
                    as.matrix(csr_dlpfc_deg_anno),
                    as.matrix(csr_dlpfc_cells_anno),
                    as.matrix(csr_dlpfc_pathway_anno))


masterTable <- rbind(as.matrix(dlpfc_pairwise2),
                     as.matrix(csr_dlpfc_pathways),
                     as.matrix(csr_dlpfc_cells),
                     as.matrix(csr_dlpfc_deg),
                     as.matrix(csr_dlpfc_ad))

write.csv(masterTable,file='~/Desktop/STG.csv',quote=F)
write.csv(masterAnno,file='~/Desktop/STGanno.csv',quote=F)
View(csr_dlpfc_pathways)
# 
# 
# mean_or<-lapply(unique(combinedScoresReducted$ModuleMethod),
#                 utilityFunctions::listify,
#                 combinedScoresReducted$aggregate,
#                 combinedScoresReducted$ModuleMethod)
# names(mean_or) <- unique(combinedScoresReducted$ModuleMethod)
# testDf<-expand.grid(names(mean_or),names(mean_or))
# testDf$pval<-mapply(function(x,y,distList){
#   return(wilcox.test(distList[[x]],
#                      distList[[y]],
#                      alternative='greater')$p.value)},
#   testDf$Var1,
#   testDf$Var2,
#   MoreArgs = list(distList=mean_or),
#   SIMPLIFY=TRUE)
# 
# testDf <- dplyr::arrange(testDf,pval)
# testDf <- testDf[1:15,]
# testDf <- dplyr::mutate(testDf,adj = p.adjust(pval,method='fdr'))
# View(testDf)
# 
# combinedScoresReducted$adjustedDegScore <- MASS::rlm(combinedScoresReducted$degScore ~ log10(combinedScoresReducted$numberOfGenes))$resid
# makeBoxPlot1(combinedScoresReducted,"adjustedDegScore","adj Number of Sig DEG enrichments")
