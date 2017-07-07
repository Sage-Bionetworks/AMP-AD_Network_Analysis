cat("\014") 


#source('loadAMPADModules.R')
source('dataPulling/pullExpressionAndPheno.R')
###get all modules from synapse now stored in allMods
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10158502")@values

###explode into list
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                           listify,
                           allMods$GeneID,
                           allMods$ModuleNameFull)
names(modulesLargeList) <- unique(allMods$ModuleNameFull)



###extract eigengenes
getEigenGenes <- function(geneList,geneExpr){
  w2 <- which(colnames(geneExpr)%in%geneList)
  if(length(w2)>10){
    geneMat <- scale(geneExpr[,w2])
    
    return(svd(geneMat)$u[,1:5])
  }else{
    return(NA)
  }
}

#get ensg genes
geneExpressionForAnalysisEnsg <- lapply(geneExpressionForAnalysis,function(x){
  y <- dplyr::select(x,dplyr::starts_with('ENSG')) %>% data.matrix
  y[is.na(y)] <- 0
  sdvec <- apply(y,2,sd)
  y <- y[,which(sdvec>0)]
  return(y)
})

#test11 <- getEigenGenes(modulesLargeList[[1]],geneExpressionForAnalysisEnsg[[1]])
test11<-utilityFunctions::outerLapply(getEigenGenes,geneExpressionForAnalysisEnsg,modulesLargeList)


adDiagnosis <- lapply(geneExpressionForAnalysis,function(x){
  return(dplyr::select(x,logitDiagnosis))
})

logisticRegressionAggregatePval <- function(x,y){
  #x <- as.numeric(x)
  y <- as.numeric(data.matrix(y))
  #print(x)
  #print(y)
  #print(str(x))
  #print(dim(y))
  logLik1 <- logLik(glm(y~x,family='binomial'))[1]
  #print(logLik1)
  logLik2 <- logLik(glm(y~1,family='binomial'))[1]
  #print(logLik2)
  lrt <- -2*logLik2+2*logLik1
  #print(lrt)
  return(pchisq(lrt,1,lower.tail=F))
}

logisticRegressionOddsRatio <- function(x,y){
  y <- as.numeric(data.matrix(y))
  glmmod <- glm(y~x,family='binomial')
  return(summary(glmmod)$coef[2,1])
}

wrapperFxn1 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn2 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      apply(eigengenes,2,logisticRegressionAggregatePval,diagnosis) %>% return
    }else{
      return(NA)
    }
  }
  lapply(eigengenes,wrapperFxn2,diagnosis) %>% return
}

wrapperFxn3 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn4 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      eigengenes = scale(eigengenes)
      apply(eigengenes,2,logisticRegressionOddsRatio,diagnosis) %>% return
    }else{
      return(NA)
    }
  }
  lapply(eigengenes,wrapperFxn4,diagnosis) %>% return
}

wrapperFxn5 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn6 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      #eigengenes = scale(eigengenes)
      cor(eigengenes,diagnosis,use = 'pairwise.complete.obs') %>% c %>% return
    }else{
      return(NA)
    }
  }
  res <- lapply(eigengenes,wrapperFxn6,diagnosis)
  names(res) <- names(eigengenes)
  return(res)
}


ModulePvalues <- mapply(wrapperFxn1,
                        adDiagnosis,
                        test11,
                        SIMPLIFY = F)

ModuleOddsRatios <- mapply(wrapperFxn3,
                           adDiagnosis,
                           test11,
                           SIMPLIFY= F)

ModuleCorrelations <- mapply(wrapperFxn5,
                             adDiagnosis,
                             test11,
                             SIMPLIFY = F)

fisherWrapper <- function(pvalues){
  fisherFxn <- function(x){
    score <- (-2*log(x)) %>% 
      sum %>%
      pchisq(df = 2*length(x),
             lower.tail=F) %>%
      return
  }
  lapply(pvalues,fisherFxn) %>%
    return
}


ModulePvalues2 <- lapply(ModulePvalues,
                         fisherWrapper)


pval1 <- unlist(ModulePvalues)
pval2 <- unlist(ModulePvalues2)

gap::qqunif(pval1)
gap::qqunif(pval2)
save.image(file = '~/backup.rda')

moduleManifest <- data.frame(ModuleNameFull = names(modulesLargeList),
                             stringsAsFactors=F)
moduleManifest <- dplyr::inner_join(moduleManifest,
                                    allMods%>%dplyr::select(Module,method,ModuleName,brainRegion,ModuleNameFull),
                                    by=c('ModuleNameFull'))
moduleManifest <- moduleManifest[which(!duplicated(moduleManifest)),]

buildSingleManifest <- function(pvalList,orList,corList,modMani){
  foo <- data.frame(pvalList,stringsAsFactors=F)
  foo <- t(foo)
  colnames(foo) <- paste0('eigengene',1:5,'pval')
  foo <- data.frame(foo,stringsAsFactors=F)
  foo$ModuleNameFull <- rownames(foo)
  print(foo[1:5,])
  
  foo_or <- data.frame(orList,stringsAsFactors = F)
  foo_or <- t(foo_or)
  colnames(foo_or) <- paste0('eigengene',1:5,'OR')
  foo_or <- data.frame(foo_or,stringsAsFactors=F)
  foo_or$ModuleNameFull <- rownames(foo_or)
  print(foo_or[1:5,])
  
  foo_cor <- data.frame(corList,stringsAsFactors = F)
  foo_cor <- t(foo_cor)
  colnames(foo_cor) <- paste0('eigengene',1:5,'cor')
  foo_cor <- data.frame(foo_cor,stringsAsFactors=F)
  foo_cor$ModuleNameFull <- rownames(foo_cor)
  print(foo_cor[1:5,])
  
  foo <- dplyr::inner_join(modMani,foo)
  print(foo[1:5,])
  foo <- dplyr::inner_join(foo,foo_or)
  print(foo[1:5,])
  foo <- dplyr::inner_join(foo,foo_cor)
  print(foo[1:5,])
  return(foo)
}

moduleManifestList <- mapply(buildSingleManifest,
                             ModulePvalues,
                             ModuleOddsRatios,
                             ModuleCorrelations,
                             MoreArgs = list(modMani=moduleManifest),SIMPLIFY=F)

moduleManifestList <- mapply(function(x,y){
  x <- dplyr::mutate(x,brainRegionAssociation = y)
  return(x)
  }, moduleManifestList,
  names(moduleManifestList),
  SIMPLIFY=FALSE)

fisherFxn <- function(x){
  x <- as.numeric(x)
  score <- (-2*log(x)) %>% 
    sum %>%
    pchisq(df = 2*length(x),
           lower.tail=F) %>%
    return
}

moduleManifestListCollapsed <- do.call(rbind,moduleManifestList)
notmissing <- which(rowSums(is.na(moduleManifestListCollapsed))==0)
moduleManifestListCollapsed <- moduleManifestListCollapsed[notmissing,]
rSynapseUtilities::makeTable(moduleManifestListCollapsed,'cross study eigengene p-values, odds ratios, pearson cor','syn2370594')

#moduleManifestCombined <- do.call(rbind,
#                                  moduleManifestList)

foo3 <- dplyr::select(moduleManifestCombined,dplyr::starts_with("eigengene"))
aggPval <- apply(foo3,1,fisherFxn)

moduleManifestCombined2 <- dplyr::mutate(moduleManifestCombined,
                                         eigengeneAggregate = aggPval)

moduleManifestCombined2 <- dplyr::arrange(moduleManifestCombined2,eigengeneAggregate)

moduleManifestCombined3 <- moduleManifestCombined2
#moduleManifestCombined3$ModuleNameFull <- gsub("kmeans","metanetwork",moduleManifestCombined3$ModuleNameFull)
#moduleManifestCombined3$method <- gsub("kmeans","metanetwork",moduleManifestCombined3$method)
#moduleManifestCombined3$ModuleName <- gsub("kmeans","metanetwork",moduleManifestCombined3$ModuleName)
rSynapseUtilities::makeTable(moduleManifestCombined3,'cross study eigengene analysis','syn2370594')










aggregated20 <- dplyr::group_by(dplyr::select(moduleManifestCombined3,ModuleNameFull,eigengeneAggregate),ModuleNameFull)
aaa1<-dplyr::summarise(aggregated20,a = sum(-log(eigengeneAggregate)))

#moduleManifestCombined3 <- gsub("kmeans","metanetwork",moduleManifestCombined2)
metanetworkMods <- dplyr::filter(moduleManifestCombined3,method=='metanetwork')





ensgrosmapDf <- dplyr::select(rosmapDf,dplyr::starts_with('ENSG')) %>% data.matrix
ensgrosmapDf[is.na(ensgrosmapDf)] <- 0
rosmapEigengenes <- lapply(modulesLargeList,getEigenGenes,ensgrosmapDf)







moduleManifest <- data.frame(ModuleNameFull = names(modulesLargeList),
                             stringsAsFactors=F)

moduleManifest <- dplyr::inner_join(moduleManifest,
                                    allMods%>%dplyr::select(Module,method,ModuleName,brainRegion,ModuleNameFull),
                                    by=c('ModuleNameFull'))

moduleSizeDistn <- sapply(modulesLargeList,length)

rosmapModuleSummaryTable <- data.frame(module=names(moduleSizeDistn),size=moduleSizeDistn,stringsAsFactors=F)
rosmapModuleSummaryTable <- dplyr::left_join(rosmapModuleSummaryTable,allMods%>%dplyr::select(ModuleName,method),c('module'='ModuleName'))
rosmapModuleSummaryTable <- rosmapModuleSummaryTable[which(!duplicated(rosmapModuleSummaryTable)),]


barplot(sort(rosmapModuleSummaryTable$size,decreasing=T),pch=15,col=rainbow(4)[as.numeric(as.factor(rosmapModuleSummaryTable$method))[order(rosmapModuleSummaryTable$size,decreasing=T)]],xlab='module',ylab='module size',main='rosmap module size distn',border=NA)
legend('topright',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))

####run association analyses in ROSMAP data for 1,2,5 eigengenes

rosmapDf <- geneExpressionForAnalysis$rosmapDLPFC

###for each module defiition produce eigengenes

getEigenGenes <- function(geneList,geneExpr){
  w2 <- which(colnames(geneExpr)%in%geneList)
  if(length(w2)>10){
    geneMat <- scale(geneExpr[,w2])
    return(svd(geneMat)$u[,1:5])
  }else{
    return(NA)
  }
}

ensgrosmapDf <- dplyr::select(rosmapDf,dplyr::starts_with('ENSG')) %>% data.matrix
ensgrosmapDf[is.na(ensgrosmapDf)] <- 0
rosmapEigengenes <- lapply(modulesLargeList,getEigenGenes,ensgrosmapDf)

diagnos <- rosmapDf$logitDiagnosis

logisticRegressionAggregatePval <- function(x,y,n=1){
  logLik1 <- logLik(glm(y~x[,n],family='binomial'))[1]
  #print(logLik1)
  logLik2 <- logLik(glm(y~1,family='binomial'))[1]
  #print(logLik2)
  lrt <- -2*logLik2+2*logLik1
  #print(lrt)
  return(pchisq(lrt,1,lower.tail=F))
}

wrapperFxn <- function(x,y,n){
  if(is.na(x)){
    return(NA)
  }else{
    return(logisticRegressionAggregatePval(x,y,n))
  }
}

foobar1 <- sapply(rosmapEigengenes,wrapperFxn,diagnos,n=1)
foobar2 <- sapply(rosmapEigengenes,wrapperFxn,diagnos,n=2)
foobar3 <- sapply(rosmapEigengenes,wrapperFxn,diagnos,n=3)
foobar4 <- sapply(rosmapEigengenes,wrapperFxn,diagnos,n=4)
foobar5 <- sapply(rosmapEigengenes,wrapperFxn,diagnos,n=5)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,rosmapEigengene1=foobar1)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,rosmapEigengene2=foobar2)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,rosmapEigengene3=foobar3)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,rosmapEigengene4=foobar4)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,rosmapEigengene5=foobar5)
rosmapModuleSummaryTable2 <- rosmapModuleSummaryTable[which(!is.na(rosmapModuleSummaryTable$rosmapEigengene1)),]
rSynapseUtilities::makeTable(rosmapModuleSummaryTable2,"Rosmap Module Summaries","syn2370594")
#a1 <- dplyr::arrange(rosmapModuleSummaryTable,rosmapEigengene1)
aa1 <- dplyr::group_by(rosmapModuleSummaryTable,method)
barplot(-log10(aa1$rosmapEigengene1),col=rainbow(4)[as.factor(aa1$method)],border=NA,xlab='modules',ylab='-log10 pval',main='principle eigengene',ylim=c(0,10))
legend('topleft',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))

barplot(-log10(aa1$rosmapEigengene2),col=rainbow(4)[as.factor(aa1$method)],border=NA,xlab='modules',ylab='-log10 pval',main='second eigengene',ylim=c(0,10))
legend('topleft',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))

barplot(-log10(aa1$rosmapEigengene3),col=rainbow(4)[as.factor(aa1$method)],border=NA,xlab='modules',ylab='-log10 pval',main='third eigengene',ylim=c(0,10))
legend('topleft',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))

barplot(-log10(aa1$rosmapEigengene4),col=rainbow(4)[as.factor(aa1$method)],border=NA,xlab='modules',ylab='-log10 pval',main='fourth eigengene',ylim=c(0,10))
legend('topleft',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))

barplot(-log10(aa1$rosmapEigengene5),col=rainbow(4)[as.factor(aa1$method)],border=NA,xlab='modules',ylab='-log10 pval',main='fifth eigengene',ylim=c(0,10))
legend('topleft',c('megena','metanetwork','speakEasy','wina'),fill=rainbow(4))


####enrichment analysis for ad genes
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList


microglialList <- GeneSets$Cell_Markers$`Zhang:Microglia`
neuralList <- GeneSets$Cell_Markers$`Zhang:Neuron`
astrocyteList <- GeneSets$Cell_Markers$`Zhang:Astrocyte`
endothelialList <- GeneSets$Cell_Markers$`Zhang:Endothelial`


masterGeneMap <- utilityFunctions::convertEnsemblToHgnc(allMods$GeneID)
allMods <- dplyr::left_join(allMods,masterGeneMap,c('GeneID'='ensembl_gene_id'))

####enrichment analysis for cell signature genes
convertEnsToHgnc <- function(x,map){
  
  return(map$external_gene_name[which(map$ensembl_gene_id%in%x)])
}
modulesLargeListHgnc <- lapply(modulesLargeList,convertEnsToHgnc,masterGeneMap)


moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined,
                                         eigengeneAggregate = aggPval)

####
adEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,adList,masterGeneMap$external_gene_name)

tdf <- data.frame(ModuleNameFull = names(unlist(adEnrich)),adGeneticEnrich = unlist(adEnrich),stringsAsFactors=F)
moduleManifestCombined2b <- dplyr::left_join(moduleManifestCombined2b,tdf)

#moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined2b,adGeneticEnrich = unlist(adEnrich))

microglialEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,microglialList,masterGeneMap$external_gene_name)
#moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined2b,microglialEnrichment = unlist(microglialEnrich))
tdf <- data.frame(ModuleNameFull = names(unlist(microglialEnrich)),microgEnrich = unlist(microglialEnrich),stringsAsFactors=F)
moduleManifestCombined2b <- dplyr::left_join(moduleManifestCombined2b,tdf)

neuralEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,neuralList,masterGeneMap$external_gene_name)
tdf <- data.frame(ModuleNameFull = names(unlist(neuralEnrich)),neurEnrich = unlist(neuralEnrich),stringsAsFactors=F)
moduleManifestCombined2b <- dplyr::left_join(moduleManifestCombined2b,tdf)


#moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined2b,neuralEnrichment = unlist(neuralEnrich))

astrocyteEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,astrocyteList,masterGeneMap$external_gene_name)
#moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined2b,astrocyteEnrichment = unlist(astrocyteEnrich))
tdf <- data.frame(ModuleNameFull = names(unlist(astrocyteEnrich)),astroEnrich = unlist(astrocyteEnrich),stringsAsFactors=F)
moduleManifestCombined2b <- dplyr::left_join(moduleManifestCombined2b,tdf)

endothelialEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,endothelialList,masterGeneMap$external_gene_name)
#moduleManifestCombined2b <- dplyr::mutate(moduleManifestCombined2b,endothelialEnrichment = unlist(endothelialEnrich))
tdf <- data.frame(ModuleNameFull = names(unlist(endothelialEnrich)),endoEnrich = unlist(endothelialEnrich),stringsAsFactors=F)
moduleManifestCombined2b <- dplyr::left_join(moduleManifestCombined2b,tdf)



View(moduleManifestCombined2b)


rSynapseUtilities::makeTable(moduleManifestCombined2b,"cross study analyses","syn2370594")
rSynapseUtilities::makeTable(allMods,"rosmap modules","syn2370594")

keepTab<-dplyr::filter(moduleManifestCombined2b,p.adjust(adGeneticEnrich,method='fdr')<=0.01)
View(keepTab)
#p.adjust(moduleManifestCombined2b$adGeneticEnrich,method='fdr')
rownames(keepTab) <- paste0(keepTab$ModuleNameFull,keepTab$brainRegionAssociation)
keepTab <- dplyr::select(keepTab,-method,-Module,-brainRegion,-brainRegionAssociation,-ModuleNameFull,-ModuleName)

keepTab2 <- keepTab<0.0001
keepTab2[which(keepTab==TRUE)] <- 1


pheatmap::pheatmap(data.matrix(keepTab2),
                   cluster_rows=F,
                   cluster_cols=F)
