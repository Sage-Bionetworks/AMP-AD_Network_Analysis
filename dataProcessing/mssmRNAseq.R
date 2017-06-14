synapseClient::synapseLogin()
#steps: 
#1) get rid of first 4 rows
#2) split cerebellum and temporal cortex samples
#3) ensure proper location for expression data and networks on Synapse done
#4) make sure output is in proper format going forward

exprDataObj <- synapseClient::synGet('syn8485027')
#geneIdObj <- synGet('syn4922926')
#sampleIdObj <- synGet('syn4922923')

exprData <- data.table::fread(exprDataObj@filePath,data.table=F)

geneKey <- exprData[,1]
exprData <- exprData[,-1]
exprData <- t(exprData)
colnames(exprData) <- geneKey
exprData <- t(exprData)


#set missing values to means
#meanNa <- rowSums(is.na(exprData))


#covariates
mssmcovObj <- synapseClient::synGet('syn6100548')
mssmcov <- data.table::fread(mssmcovObj@filePath,data.table=F)
library(dplyr)
fpIds <- dplyr::filter(mssmcov,BrodmannArea=='BM10')%>%
  dplyr::select(sampleIdentifier)
stgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM22')%>%
  dplyr::select(sampleIdentifier)
phgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM36')%>%
  dplyr::select(sampleIdentifier)
IFGIds <- dplyr::filter(mssmcov,BrodmannArea=='BM44')%>%
  dplyr::select(sampleIdentifier)

exprDataFP <- exprData[,colnames(exprData)%in%fpIds$sampleIdentifier]
exprDataSTG <- exprData[,colnames(exprData)%in%stgIds$sampleIdentifier]
exprDataPHG <- exprData[,colnames(exprData)%in%phgIds$sampleIdentifier]
exprDataIFG <- exprData[,colnames(exprData)%in%IFGIds$sampleIdentifier]

pushToSynapseFxn <- function(exprMat,fileName,parentId,comment){
  write.csv(exprMat,file=fileName,quote=F)
  foo = synapseClient::File(fileName,parentId=parentId,versionComment = comment)
  bar2 = synapseClient::synGet('syn5898491',downloadFile=F)
  anno = synapseClient::synGetAnnotations(bar2)
  synapseClient::synSetAnnotations(foo) = as.list(anno)
  permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                                 ref = 'branch',
                                 refName = 'module-comparisons',
                                 repositoryPath = 'mssmRNAseq.R')
  foo = synapseClient::synStore(foo,
                 used = as.list(c('syn8485027')),
                 executed = as.list(c(permLink)),
                 activityName = 'Format MSSM RNAseq Data',
                 activityDescription = 'Push MSSM RNAseq data into format for network pipeline')
}

pushToSynapseFxn(exprDataFP,'MSSMRNAseq_FP.csv',"syn8281257","MSSM RNAseq FP expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia with no dx adj for May 2017 data freeze")
pushToSynapseFxn(exprDataSTG,'MSSMRNAseq_STG.csv',"syn8281286","MSSM RNAseq STG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia with no dx adj for May 2017 data freeze")
pushToSynapseFxn(exprDataPHG,'MSSMRNAseq_PHG.csv',"syn8281279","MSSM RNAseq PHG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia with no dx adj for May 2017 data freeze")
pushToSynapseFxn(exprDataIFG,'MSSMRNAseq_IFG.csv',"syn8281272","MSSM RNAseq IFG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia with no dx adj for May 2017 data freeze")
