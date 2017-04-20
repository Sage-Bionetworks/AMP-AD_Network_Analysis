synapseClient::synapseLogin()
rosmapClinicalObj <- synapseClient::synGet('syn3191087')
rosmapUncensoredAgesObj <- synapseClient::synGet('syn7116000')
rosmapIdMapObj <- synapseClient::synGet('syn3382527')
rosmapCogDecline1Obj <- synapseClient::synGet('syn6182375')
rosmapCogDecline2Obj <- synapseClient::synGet('syn6182376')
rosmapMetanetworkObj <- synapseClient::synGet('syn8268669')

####thanneer's code for fixing ids
#https://github.com/th1vairam/Brain_Reg_Net/blob/ad12b544d6d6c23be2bf94eed53a6b4f75d154d1/code/Rmd/ROSMAP_REPROCESSED.Rmd

####read in clinical file
rosmapClinical <- data.table::fread(synapseClient::getFileLocation(rosmapClinicalObj),
                                    data.table=F)
rosmapClinical <- dplyr::select(rosmapClinical,-V1,-age_death,-age_first_ad_dx,-age_at_visit_max)
#View(rosmapClinical)


####read in uncensored age file
rosmapUncensoredAges <- data.table::fread(synapseClient::getFileLocation(rosmapUncensoredAgesObj),
                                          data.table=F)
#View(rosmapUncensoredAges)

####merge uncensored ages with clinical file
rosmapClinical <- dplyr::left_join(rosmapClinical,rosmapUncensoredAges,by='projid')
#View(rosmapClinical)

###cognitive decline
rosmapCogDec1 <- data.table::fread(synapseClient::getFileLocation(rosmapCogDecline1Obj),
                                   data.table=F)
rosmapCogDec1 <- dplyr::select(rosmapCogDec1,-Sample)

rosmapCogDec2 <- data.table::fread(synapseClient::getFileLocation(rosmapCogDecline2Obj),
                                   data.table=F)

mungeIds <- function(x){
  #ros or map
  isROS<-grep('ROS',x)
  isMAP<-grep('MAP',x)
  if(length(isROS)>0){
    rosId<-strsplit(x,'ROS')[[1]][2]
    return(as.numeric(rosId))
  }else if(length(isMAP)>0){
    mapId <- strsplit(x,'MAP')[[1]][2]
    return(as.numeric(mapId))
  }else{
    return(as.numeric(x))
  }
}
projids2 <- sapply(rosmapCogDec2$CollaboratorParticipantId,mungeIds)
rosmapCogDec2$ProjectID <- projids2
rosmapCogDec2 <- dplyr::select(rosmapCogDec2,ProjectID,cogn_global_slope)
rosmapCogDec <- rbind(rosmapCogDec1,rosmapCogDec2)

#View(rosmapCogDec1)
rosmapClinical <- dplyr::left_join(rosmapClinical,
                                   rosmapCogDec,
                                   by=c("projid"="ProjectID"))
rosmapClinical$apoe_genotype <- as.factor(rosmapClinical$apoe_genotype)
rosmapClinical$ceradsc <- as.factor(rosmapClinical$ceradsc)
rosmapClinical$braaksc <- as.factor(rosmapClinical$braaksc)
rosmapClinical$cogdx <- as.factor(rosmapClinical$cogdx)
####simple tests
lmObj <- lm(cogn_global_slope ~ age_death+ceradsc+braaksc+cogdx+apoe_genotype,data = rosmapClinical)
summary(lmObj)
eff1<-effects::allEffects(lmObj)

library(dplyr)
KEY <- read.csv(synapseClient::getFileLocation(rosmapIdMapObj))%>%
  dplyr::filter(mrna_data == 1) %>%
  dplyr::select(projid, mrna_id) %>%
  tidyr::separate(mrna_id, c('a','b','batch'), sep = '_') %>%
  tidyr::unite(Sampleid, a, b) %>%
  dplyr::select(-batch) %>%
  unique

expressionDataObj <- synapseClient::synGet('syn8456719')
expressionData <- data.table::fread(synapseClient::getFileLocation(expressionDataObj),data.table=F)
rownames(expressionData) <- expressionData$ensembl_gene_id
expressionData <- dplyr::select(expressionData,-ensembl_gene_id)
expressionData <- t(expressionData)
expressionData <- data.frame(expressionData,stringsAsFactors=F)
expressionData$aSampleId <- rownames(expressionData)

expressionData <- dplyr::left_join(expressionData,KEY,by=c('aSampleId'='Sampleid'))
combinedData <- dplyr::left_join(expressionData,rosmapClinical,by='projid')



unlist(test11)
names(baz4[[1]])

#lmObj <- lm(cogn_global_slope ~ age_death+ceradsc+braaksc+cogdx+apoe_genotype+ENSG00000182621,data = combinedData)
#summary(lmObj)
load(synapseClient::getFileLocation(rosmapMetanetworkObj))

rosmap <- igraph::graph_from_adjacency_matrix(bicNetworks$network,
                                              mode='undirected')
degree <- igraph::degree(rosmap)
degree2 <- data.frame(ensembl_gene_id=names(degree),degree=degree,stringsAsFactors=F)
degree2 <- dplyr::arrange(degree2,desc(degree))
degree2[1:5,]

foobar <- model.matrix(~1 + age_death+ceradsc+braaksc+cogdx+apoe_genotype,data=combinedData)
foobar2 <- model.matrix(~1,data=combinedData)
a1 <-which(is.na(combinedData$cogn_global_slope))
foobar <- foobar[-a1,]
foobar2 <- foobar2[-a1,]
combinedData <- combinedData[-a1,]


expressionHarmonized <- dplyr::select(combinedData,starts_with("ENSG"))

var2b <- utilityFunctions::varCompWrapperFunction(outcome=combinedData$cogn_global_slope,
                                                  fixed = foobar2,
                                                  randomFeatures=expressionHarmonized)
var1b <- utilityFunctions::varCompWrapperFunction(outcome=combinedData$cogn_global_slope,
                                                  fixed = foobar,
                                                  randomFeatures = expressionHarmonized)

barplot(c(var2b$percentVarianceExplained,var1b$percentVarianceExplained),
        col=c('red','blue'),
        names.arg=c('~1','~1+age+cerad+braak+cogdx+apoe'),
        xlab='model',
        ylab='% variation explained',
        main='% variation explained for rate of cog decline by expression',
        ylim=c(0,0.5))
text(c(.7,1.9),c(0.48,0.1),c('p<1e-16','p=0.03'))
#utilityFunctions::varCompWrapperFunction()
load('rosmapMods.rda')
#module109genes <- synapseClient::synTableQuery("SELECT geneName FROM syn5321231 where speakeasyModule=109")@values
#baz4[[1]]$module109 <- unique(module109genes$geneName)
getVariance <- function(moduleGenes,combinedData1,fixedData1){
  foo1 <- utilityFunctions::varCompWrapperFunction(outcome=combinedData1$cogn_global_slope,
                                                   fixed = fixedData1,
                                                   randomFeatures = combinedData1[,which(colnames(combinedData1)%in%moduleGenes)])
  return(unlist(foo1))
}
test11 <- sapply(abcd[[3]],getVariance,combinedData1=combinedData,fixedData1=foobar2)
test11 <- t(test11)
test11 <- data.frame(test11,stringsAsFactors=F)
test11$module <- rownames(test11)
test11 <- dplyr::arrange(test11,pvalue)
test11$fdr <- p.adjust(test11$pvalue,method='fdr')
View(test11)

multiModelKernels <- function(geneList,x){
  return(cor(t(x[,geneList]),use = 'pairwise.complete.obs'))
}
getKernels <- lapply(baz4[[1]],multiModelKernels,combinedData)

fullMod <- varComp::varComp(combinedData$cogn_global_slope ~ foobar,varcov=cor(t(scale(expressionHarmonized))))
fullMod <- varComp::varComp(combinedData$cogn_global_slope ~ foobar-1,varcov=getKernels)

test111 <- test11[-which(test11$module=="module109"),]
View(test111)
par(mar=c(6,4,4,4))
barplot(test111$percentVarianceExplained*100,
        names.arg = test111$module,
        col = test111$module,
        xlab='module',
        ylab='% variation explained',
        main='rate of cognitive decline',
        las=2,
        ylim=c(0,6.2))
text("***",x = (c(1,2,3,4,5)-.5+(1:5)*.2),y=test111$percentVarianceExplained[1:5]*100+.2)
legend('topright',c("*** FDR<=0.05"),cex=.75,box.lty=0)

gap::qqunif(test11[,'pvalue'])

lmFun <- function(x,y,z){
  lmObj <- lm(y~x+z-1)
  return(summary(lmObj)$coef[1,4])
}

#expressionHarmonized <- expressionHarmonized[-a1,]
foo2 <- apply(expressionHarmonized,
              2,
              lmFun,
              y=as.matrix(combinedData$cogn_global_slope),
              z=foobar)
#nullModel <- utilityFunctions::varCompWrapperFunction(outcome=scale(combinedData$cogn_global_slope),
#                                                     features=data.matrix(foobar))

tabl23 <- data.frame(pval=foo2,ensembl_gene_id=names(foo2),stringsAsFactors = F)
#merge tabl23 and degree2
degree3 <- dplyr::left_join(degree2,tabl23,by='ensembl_gene_id')
degree3 <- dplyr::arrange(degree3,desc(degree))

###lasso
lassoFxn <- function(foobar,expressionHarmonized,aleph=1,nfolds1=10){
  masterMatrix <- cbind(data.matrix(foobar),data.matrix(expressionHarmonized))
  set.seed(1)
  lasso <- glmnet::cv.glmnet(y=scale(combinedData$cogn_global_slope),
                             x=masterMatrix,
                             penalty.factor=c(rep(0,ncol(foobar)),rep(1,ncol(expressionHarmonized))),
                             alpha=aleph,
                             nfolds=nfolds1)
  model <- list()
  model$min_error <- min(lasso$cvm)
  model$path_min <- which.min(lasso$cvm)
  model$coef <- which(lasso$glmnet.fit$beta[,model$path_min]!=0)
  model$lasso <- lasso
  return(model)
}
lassoAlt1 <- lassoFxn(foobar,expressionHarmonized)
lassoAlt2 <- lassoFxn(foobar,dplyr::select(expressionHarmonized,one_of(degree3$ensembl_gene_id[1:1500])),aleph=1);c(lassoAlt2$min_error,lassoAlt2$path_min)



collateGraphStatistics <- function(graph){
  model <- list()
  cat('Computing Degree...\n')
  model$degree <- igraph::degree(graph)
  #model$alpha_centrality <- igraph::alpha_centrality(graph,tol=1e-14)
  cat('Computing Authority Score...\n')
  model$authority_score <- igraph::authority_score(graph)$vector
  
  cat('Computing Closeness...\n')
  model$closeness <- igraph::closeness(graph)
  
  cat('Computing Eccentricity...\n')
  model$eccentricity <- igraph::eccentricity(graph)
  
  cat('Computing Eigenvector Centrality...\n')
  model$eigen_centrality <- igraph::eigen_centrality(graph)$vector
  
  cat('Computing Betweeness Centrality...\n')
  model$centr_betw <- igraph::betweenness(graph)
  
  cat('Computing Transitivity...\n')
  model$transitivity <- igraph::transitivity(graph,
                                             type='undirected')
  
  return(model)
}

rosmapNetworkStats <- collateGraphStatistics(rosmap)
rosmapNetworkStats2 <- data.frame(rosmapNetworkStats,stringsAsFactors=F)
rosmapNetworkStats2$ensembl_gene_id <- names(degree)
rosmapNetworkStats2 <- dplyr::arrange(rosmapNetworkStats2,desc(authority_score))
tabl34 <- dplyr::left_join(rosmapNetworkStats2,tabl23,by='ensembl_gene_id')
tabl34 <- dplyr::arrange(tabl34,desc(closeness))
gap::qqunif(tabl34$pval,xlim=c(0,5),ylim=c(0,6))
par(new=T)
gap::qqunif(tabl34$pval[1:30],col='red',xlim=c(0,5),ylim=c(0,6))

lassoAlt1 <- lassoFxn(foobar,expressionHarmonized)
lassoAlt2 <- lassoFxn(foobar,dplyr::select(expressionHarmonized,one_of(tabl34$ensembl_gene_id[1:30])),aleph=1,nfolds=510);c(lassoAlt2$min_error,lassoAlt2$path_min)


geneMap <- utilityFunctions::convertEnsemblToHgnc(tabl34$ensembl_gene_id)
tabl56 <- dplyr::left_join(tabl34,geneMap,by='ensembl_gene_id')
View(tabl56)
y <- -log10(tabl56$pval)
summary(lm(y~degree+authority_score+closeness+eccentricity+eigen_centrality+centr_betw+transitivity,data=tabl56))
betas <- summary(a11)$coef[,1]
tabl56 <- dplyr::mutate(tabl56,networkFeature=betas[1] + betas[2]*degree + betas[3]*closeness + betas[4]*eccentricity + betas[5]*eigen_centrality + betas[6]*centr_betw)
tabl56 <- dplyr::arrange(tabl56,desc(networkFeature))
tabl56 <- dplyr::arrange(tabl56,(pval))
cat(tabl56$external_gene_name[1:500],file='~/Desktop/networkz.csv',sep='\n')

