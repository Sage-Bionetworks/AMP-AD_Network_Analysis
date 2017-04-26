synapseClient::synapseLogin()

#get synIds for gene expression variables
geneExpressionDataManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( dataSubType = 'residualGeneExpForNetAnlz' ) AND ( normalizationType = 'CQN' ) )")

#get synIds for covariates
covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

#geneExpressionDataObj <- sapply(geneExpressionDataManifest@values$id,synapseClient::synGet)
#covariateManifestObj <- sapply(covariateManifest@values$id,synapseClient::synGet)

#load expression data into R
geneExpressionList <- rSynapseUtilities::loadDelimIntoList(geneExpressionDataManifest@values$id)

#load covariate data into R
covariateList <- rSynapseUtilities::loadDelimIntoList(covariateManifest@values$id)


#split mayo into two data-frames
geneExpressionForAnalysis <- list()
geneExpressionForAnalysis$mayoTCX <- dplyr::select(geneExpressionList$syn8466826,
                         dplyr::ends_with('TCX'))
rownames(geneExpressionForAnalysis$mayoTCX) <- geneExpressionList$syn8466826$ensembl_gene_id

geneExpressionForAnalysis$mayoCER <- dplyr::select(geneExpressionList$syn8466826,
                         dplyr::ends_with('CER'))
rownames(geneExpressionForAnalysis$mayoCER) <- geneExpressionList$syn8466826$ensembl_gene_id
library(dplyr)
#rosmap
geneExpressionForAnalysis$rosmapDLPFC <- geneExpressionList$syn8456719 %>% dplyr::select(-ensembl_gene_id)
rownames(geneExpressionForAnalysis$rosmapDLPFC) <- geneExpressionList$syn8456719$ensembl_gene_id
#split mssm into 4 data-frames
mssmcovObj <- synapseClient::synGet('syn6100548')
mssmcov <- data.table::fread(mssmcovObj@filePath,data.table=F)
fpIds <- dplyr::filter(mssmcov,BrodmannArea=='BM10')%>%
  dplyr::select(sampleIdentifier)
stgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM22')%>%
  dplyr::select(sampleIdentifier)
phgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM36')%>%
  dplyr::select(sampleIdentifier)
IFGIds <- dplyr::filter(mssmcov,BrodmannArea=='BM44')%>%
  dplyr::select(sampleIdentifier)

geneExpressionForAnalysis$msbbFP <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%fpIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbSTG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%stgIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbPHG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%phgIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbIFG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%IFGIds$sampleIdentifier)]
rownames(geneExpressionForAnalysis$msbbFP) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbSTG) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbPHG) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbIFG) <- geneExpressionList$syn8485027$ensembl_gene_id

####transpose all matrices
geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,t)

####make data frames
geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,data.frame,stringsAsFunctions=F)

####add sample id as first column
addSampleId <- function(x){
  x <- dplyr::mutate(x,aSampleId=rownames(x))
  return(x)
}
geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,addSampleId)

 
#add diagnosis to each expression data frame with left join
#mayo
geneExpressionForAnalysis$mayoTCX <- dplyr::left_join(geneExpressionForAnalysis$mayoTCX,
                                                      covariateList$syn8466814%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                      c('aSampleId'='SampleID'))
w1<-which(colnames(geneExpressionForAnalysis$mayoTCX)%in%c('aSampleId','BrainRegion.Diagnosis'))
otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$mayoTCX),w1)
geneExpressionForAnalysis$mayoTCX <- geneExpressionForAnalysis$mayoTCX[,c(w1,otherCol)]
colnames(geneExpressionForAnalysis$mayoTCX)[c(1,2)] <- c('SampleID','Diagnosis')

logitDiag <- sapply(geneExpressionForAnalysis$mayoTCX$Diagnosis,function(x){
  if(x=='TCX.AD'){
    return(1)
  }else if (x =='TCX.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
mayoTCX <- dplyr::mutate(geneExpressionForAnalysis$mayoTCX,
                                                      logitDiagnosis = logitDiag)

geneExpressionForAnalysis$mayoTCX <- mayoTCX




geneExpressionForAnalysis$mayoCER <- dplyr::left_join(geneExpressionForAnalysis$mayoCER,
                                                      covariateList$syn8466814%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                      c('aSampleId'='SampleID'))
w1<-which(colnames(geneExpressionForAnalysis$mayoCER)%in%c('aSampleId','BrainRegion.Diagnosis'))
otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$mayoCER),w1)
geneExpressionForAnalysis$mayoCER <- geneExpressionForAnalysis$mayoCER[,c(w1,otherCol)]
colnames(geneExpressionForAnalysis$mayoCER)[c(1,2)] <- c('SampleID','Diagnosis')

logitDiag <- sapply(geneExpressionForAnalysis$mayoCER$Diagnosis,function(x){
  if(x=='CER.AD'){
    return(1)
  }else if (x =='CER.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
mayoCER <- dplyr::mutate(geneExpressionForAnalysis$mayoCER,
                         logitDiagnosis = logitDiag)
geneExpressionForAnalysis$mayoCER <- mayoCER


###rosmap
geneExpressionForAnalysis$rosmapDLPFC <- dplyr::left_join(geneExpressionForAnalysis$rosmapDLPFC,
                                                        covariateList$syn8456631%>%dplyr::select(SampleID,Diagnosis),
                                                        c('aSampleId'='SampleID'))
w1<-which(colnames(geneExpressionForAnalysis$rosmapDLPFC)%in%c('aSampleId','Diagnosis'))
otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$rosmapDLPFC),w1)
geneExpressionForAnalysis$rosmapDLPFC <- geneExpressionForAnalysis$rosmapDLPFC[,c(w1,otherCol)]
colnames(geneExpressionForAnalysis$rosmapDLPFC)[c(1,2)] <- c('SampleID','Diagnosis')

logitDiag <- sapply(geneExpressionForAnalysis$rosmapDLPFC$Diagnosis,function(x){
  if(x=='AD'){
    return(1)
  }else if (x =='CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
geneExpressionForAnalysis$rosmapDLPFC <- dplyr::mutate(geneExpressionForAnalysis$rosmapDLPFC,
                                                       logitDiagnosis = logitDiag)


####mssm

geneExpressionForAnalysis$msbbFP <- dplyr::left_join(geneExpressionForAnalysis$msbbFP,
                                                     covariateList$syn8484996%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                     c('aSampleId'='SampleID'))


geneExpressionForAnalysis$msbbSTG <- dplyr::left_join(geneExpressionForAnalysis$msbbSTG,
                                                     covariateList$syn8484996%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                     c('aSampleId'='SampleID'))


geneExpressionForAnalysis$msbbPHG <- dplyr::left_join(geneExpressionForAnalysis$msbbPHG,
                                                     covariateList$syn8484996%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                     c('aSampleId'='SampleID'))


geneExpressionForAnalysis$msbbIFG <- dplyr::left_join(geneExpressionForAnalysis$msbbIFG,
                                                     covariateList$syn8484996%>%dplyr::select(SampleID,BrainRegion.Diagnosis),
                                                     c('aSampleId'='SampleID'))


logitDiag <- sapply(geneExpressionForAnalysis$msbbFP$BrainRegion.Diagnosis,function(x){
  if(x=='FP.AD'){
    return(1)
  }else if (x =='FP.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
geneExpressionForAnalysis$msbbFP <- dplyr::mutate(geneExpressionForAnalysis$msbbFP,
                                                       logitDiagnosis = logitDiag)


logitDiag <- sapply(geneExpressionForAnalysis$msbbSTG$BrainRegion.Diagnosis,function(x){
  if(x=='STG.AD'){
    return(1)
  }else if (x =='STG.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
geneExpressionForAnalysis$msbbSTG <- dplyr::mutate(geneExpressionForAnalysis$msbbSTG,
                                                  logitDiagnosis = logitDiag)

logitDiag <- sapply(geneExpressionForAnalysis$msbbPHG$BrainRegion.Diagnosis,function(x){
  if(x=='PHG.AD'){
    return(1)
  }else if (x =='PHG.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
geneExpressionForAnalysis$msbbPHG <- dplyr::mutate(geneExpressionForAnalysis$msbbPHG,
                                                   logitDiagnosis = logitDiag)

logitDiag <- sapply(geneExpressionForAnalysis$msbbIFG$BrainRegion.Diagnosis,function(x){
  if(x=='IFG.AD'){
    return(1)
  }else if (x =='IFG.CONTROL'){
    return(0)
  }else{
    return(NA)
  }
})
geneExpressionForAnalysis$msbbIFG <- dplyr::mutate(geneExpressionForAnalysis$msbbIFG,
                                                   logitDiagnosis = logitDiag)



###get all modules from synapse

moduleTable <- synapseClient::synTableQuery('SELECT * FROM syn9705028')@values

###explode into list
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
modulesLargeList <- lapply(unique(moduleTable$ModuleName),
                           listify,
                           moduleTable$GeneID,
                           moduleTable$ModuleName)
names(modulesLargeList) <- unique(moduleTable$ModuleName)


moduleSizeDistn <- sapply(modulesLargeList,length)

rosmapModuleSummaryTable <- data.frame(module=names(moduleSizeDistn),size=moduleSizeDistn,stringsAsFactors=F)
rosmapModuleSummaryTable <- dplyr::left_join(rosmapModuleSummaryTable,moduleTable%>%dplyr::select(ModuleName,method),c('module'='ModuleName'))
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


masterGeneMap <- utilityFunctions::convertEnsemblToHgnc(moduleTable$GeneID)
moduleTable <- dplyr::left_join(moduleTable,masterGeneMap,c('GeneID'='ensembl_gene_id'))

####enrichment analysis for cell signature genes
convertEnsToHgnc <- function(x,map){
  
  return(map$external_gene_name[which(map$ensembl_gene_id%in%x)])
}
modulesLargeListHgnc <- lapply(modulesLargeList,convertEnsToHgnc,masterGeneMap)

####
adEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,adList,masterGeneMap$external_gene_name)
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,adGeneticEnrich = unlist(adEnrich))

microglialEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,microglialList,masterGeneMap$external_gene_name)
neuralEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,neuralList,masterGeneMap$external_gene_name)
astrocyteEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,astrocyteList,masterGeneMap$external_gene_name)
endothelialEnrich <- lapply(modulesLargeListHgnc,utilityFunctions::fisherWrapperPval,endothelialList,masterGeneMap$external_gene_name)

rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,microglialEnrichment = unlist(microglialEnrich))
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,neuralEnrichment = unlist(neuralEnrich))
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,astrocyteEnrichment = unlist(astrocyteEnrich))
rosmapModuleSummaryTable <- dplyr::mutate(rosmapModuleSummaryTable,endothelialEnrichment = unlist(endothelialEnrich))
View(rosmapModuleSummaryTable)


rSynapseUtilities::makeTable(rosmapModuleSummaryTable,"rosmap module summaries","syn2370594")
rSynapseUtilities::makeTable(moduleTable,"rosmap modules","syn2370594")
keepTab<-dplyr::filter(rosmapModuleSummaryTable,p.adjust(adGeneticEnrich,method='fdr')<=0.05)
#p.adjust(rosmapModuleSummaryTable$adGeneticEnrich,method='fdr')
rownames(keepTab) <- keepTab$module
keepTab <- dplyr::select(keepTab,-size,-method,-module)

keepTab2 <- keepTab<0.0001
keepTab2[which(keepTab)] <- 1
keepTab2[which(!keepTab)] <- 0

pheatmap::pheatmap((keepTab2),
                   cluster_rows=F,
                   cluster_cols=F)
