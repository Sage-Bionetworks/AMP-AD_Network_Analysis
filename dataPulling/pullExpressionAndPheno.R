synapseClient::synapseLogin()

#get synIds for gene expression variables
geneExpressionDataManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( dataSubType = 'residualGeneExpForNetAnlz' ) AND ( normalizationType = 'CQN' ) AND ( center <> 'ALL'))")

#get synIds for covariates
covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

#geneExpressionDataObj <- sapply(geneExpressionDataManifest@values$id,synapseClient::synGet)
#covariateManifestObj <- sapply(covariateManifest@values$id,synapseClient::synGet)

#load expression data into R
geneExpressionList <- rSynapseUtilities::loadDelimIntoList(geneExpressionDataManifest@values$id)

#load covariate data into R
covariateList <- rSynapseUtilities::loadDelimIntoList(covariateManifest@values$id)
mssmcovObj <- synapseClient::synGet('syn6100548')

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
geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,data.frame,stringsAsFactors=FALSE)

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
#w1<-which(colnames(geneExpressionForAnalysis$mayoTCX)%in%c('aSampleId','BrainRegion.Diagnosis'))
#otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$mayoTCX),w1)
#geneExpressionForAnalysis$mayoTCX <- geneExpressionForAnalysis$mayoTCX[,c(w1,otherCol)]
#colnames(geneExpressionForAnalysis$mayoTCX)[c(1,2)] <- c('SampleID','Diagnosis')

logitDiag <- sapply(geneExpressionForAnalysis$mayoTCX$BrainRegion.Diagnosis,function(x){
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
#w1<-which(colnames(geneExpressionForAnalysis$mayoCER)%in%c('aSampleId','BrainRegion.Diagnosis'))
#otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$mayoCER),w1)
#geneExpressionForAnalysis$mayoCER <- geneExpressionForAnalysis$mayoCER[,c(w1,otherCol)]
#colnames(geneExpressionForAnalysis$mayoCER)[c(1,2)] <- c('SampleID','Diagnosis')

logitDiag <- sapply(geneExpressionForAnalysis$mayoCER$BrainRegion.Diagnosis,function(x){
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
#w1<-which(colnames(geneExpressionForAnalysis$rosmapDLPFC)%in%c('aSampleId','Diagnosis'))
#otherCol <- setdiff(1:ncol(geneExpressionForAnalysis$rosmapDLPFC),w1)
#geneExpressionForAnalysis$rosmapDLPFC <- geneExpressionForAnalysis$rosmapDLPFC[,c(w1,otherCol)]
#colnames(geneExpressionForAnalysis$rosmapDLPFC)[c(1,2)] <- c('SampleID','Diagnosis')

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


