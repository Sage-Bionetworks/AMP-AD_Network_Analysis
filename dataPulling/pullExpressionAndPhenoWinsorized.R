synapseClient::synapseLogin()
library(dplyr)
#get synIds for gene expression variables
geneExpressionDataManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 where dataType = 'mRNA' and columnScaled = 'TRUE'")

#get synIds for covariates
covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

#geneExpressionDataObj <- sapply(geneExpressionDataManifest@values$id,synapseClient::synGet)
#covariateManifestObj <- sapply(covariateManifest@values$id,synapseClient::synGet)

#load expression data into R
geneExpressionList <- rSynapseUtilities::loadDelimIntoList(geneExpressionDataManifest@values$id)

#load covariate data into R
covariateList <- rSynapseUtilities::loadDelimIntoList(covariateManifest@values$id)
covariateList$syn8484996 <- synapseClient::synGet('syn8484996',version = 8) %>%
  synapseClient::getFileLocation() %>%
  data.table::fread(data.table=F)

covariateList$syn8456631 <- synapseClient::synGet('syn8456631',version = 15) %>%
  synapseClient::getFileLocation() %>%
  data.table::fread(data.table=F)

covariateList$syn8466814 <- synapseClient::synGet('syn8466814',version = 12) %>%
  synapseClient::getFileLocation() %>%
  data.table::fread(data.table=F)

#split mayo into two data-frames
geneExpressionForAnalysis <- list()
geneExpressionForAnalysis$mayoTCX <- geneExpressionList$syn8303274
colnames(geneExpressionForAnalysis$mayoTCX)[1] <- 'aSampleId'


geneExpressionForAnalysis$mayoCER <- geneExpressionList$syn8303281
colnames(geneExpressionForAnalysis$mayoCER)[1] <- 'aSampleId'
library(dplyr)
#rosmap
geneExpressionForAnalysis$rosmapDLPFC <- geneExpressionList$syn8303260
colnames(geneExpressionForAnalysis$rosmapDLPFC)[1] <- 'aSampleId'
#msbb
geneExpressionForAnalysis$msbbFP <- geneExpressionList$syn8303298
colnames(geneExpressionForAnalysis$msbbFP)[1] <- 'aSampleId'

geneExpressionForAnalysis$msbbSTG <- geneExpressionList$syn8303308
colnames(geneExpressionForAnalysis$msbbSTG)[1] <- 'aSampleId'

geneExpressionForAnalysis$msbbPHG <- geneExpressionList$syn8303314
colnames(geneExpressionForAnalysis$msbbPHG)[1] <- 'aSampleId'

geneExpressionForAnalysis$msbbIFG <- geneExpressionList$syn8303338
colnames(geneExpressionForAnalysis$msbbIFG)[1] <- 'aSampleId'


####transpose all matrices
#geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,t)

####make data frames
geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,data.frame,stringsAsFactors=FALSE)

####add sample id as first column
#addSampleId <- function(x){
#  x <- dplyr::mutate(x,aSampleId=rownames(x))
#  return(x)
#}
#geneExpressionForAnalysis <- lapply(geneExpressionForAnalysis,addSampleId)

 
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


