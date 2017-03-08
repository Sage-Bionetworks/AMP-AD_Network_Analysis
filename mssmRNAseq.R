library(synapseClient)
library(githubr)
synapseLogin()



#steps: 
#1) get rid of first 4 rows
#2) split cerebellum and temporal cortex samples
#3) ensure proper location for expression data and networks on Synapse done
#4) make sure output is in proper format going forward

exprDataObj <- synGet('syn8049659')
#geneIdObj <- synGet('syn4922926')
#sampleIdObj <- synGet('syn4922923')



populateFolders <- function(synId){
  objectsList <- list()
  objectsList$exprFolder <- synapseClient::Folder(name="Expression Data",
                                      parentId=synId)
  objectsList$exprFolder <- synapseClient::synStore(objectsList$exprFolder)
  objectsList$exprNetworks <- synapseClient::Folder(name="Expression Networks",
                                        parentId=synId)
  objectsList$exprNetworks <- synapseClient::synStore(objectsList$exprNetworks)
  objectsList$metanetwork <- synapseClient::Folder(name='metanetwork',parentId=objectsList$exprNetworks@properties$id)
  objectsList$metanetwork <- synapseClient::synStore(objectsList$metanetwork)
  objectsList$speakeasy <- synapseClient::Folder(name='speakeasy',parentId=objectsList$exprNetworks@properties$id)
  objectsList$speakeasy <- synapseClient::synStore(objectsList$speakeasy)
  objectsList$trena <- synapseClient::Folder(name='trena',parentId=objectsList$exprNetworks@properties$id)
  objectsList$trena <- synapseClient::synStore(objectsList$trena)
  objectsList$megena <- synapseClient::Folder(name='megena',parentId=objectsList$exprNetworks@properties$id)
  objectsList$megena <- synapseClient::synStore(objectsList$megena)
  objectsList$wgcna <- synapseClient::Folder(name='wgcna',parentId=objectsList$exprNetworks@properties$id)
  objectsList$wgcna <- synapseClient::synStore(objectsList$wgcna)
  return(objectsList)
}


MSSMFPId <- 'syn7980638'
MSSMITGId <- 'syn7980629'
MSSMPHGId <- 'syn7980634'
MSSMSTGId <- 'syn7980616'
MSSMFPManifest <- populateFolders(MSSMFPId)
MSSMITGManifest <- populateFolders(MSSMITGId)
MSSMPHGManifest <- populateFolders(MSSMPHGId)
MSSMSTGManifest <- populateFolders(MSSMSTGId)

library(data.table)
exprData <- fread(exprDataObj@filePath,data.table=F)
exprData <- exprData[-c(1:4),]
whichDup <- which(duplicated(exprData$ensembl_gene_id))
exprData <- exprData[-whichDup,]


geneKey <- exprData[,1:3]
exprData <- exprData[,-c(1:3)]
exprData <- t(exprData)
colnames(exprData) <- geneKey$ensembl_gene_id
exprData <- t(exprData)

#covariates
mssmcovObj <- synGet('syn6100548')
mssmcov <- fread(mssmcovObj@filePath,data.table=F)
library(dplyr)
fpIds <- dplyr::filter(mssmcov,BrodmannArea=='BM10')%>%
  dplyr::select(sampleIdentifier)
stgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM22')%>%
  dplyr::select(sampleIdentifier)
phgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM36')%>%
  dplyr::select(sampleIdentifier)
itgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM44')%>%
  dplyr::select(sampleIdentifier)



exprDataFP <- exprData[,colnames(exprData)%in%fpIds$sampleIdentifier]
exprDataSTG <- exprData[,colnames(exprData)%in%stgIds$sampleIdentifier]
exprDataPHG <- exprData[,colnames(exprData)%in%phgIds$sampleIdentifier]
exprDataITG <- exprData[,colnames(exprData)%in%itgIds$sampleIdentifier]

pushToSynapseFxn <- function(exprMat,fileName,parentId,comment){
  write.csv(exprMat,file=fileName,quote=F)
  foo = File(fileName,parentId=parentId,versionComment = comment)
  bar2 = synGet('syn5898491',downloadFile=F)
  anno = synGetAnnotations(bar2)
  synSetAnnotations(foo) = as.list(anno)
  permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                                 ref = 'branch',
                                 refName = 'mssm',
                                 repositoryPath = 'mssmRNAseq.R')
  foo = synStore(foo,
                 used = as.list(c('syn8049659')),
                 executed = as.list(c(permLink)),
                 activityName = 'Format MSSM RNAseq Data',
                 activityDescription = 'Push MSSM RNAseq data into format for network pipeline')
}

pushToSynapseFxn(exprDataFP,'MSSMRNAseq_FP.csv',MSSMFPManifest$exprFolder@properties$id,"MSSM RNAseq FP expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows")
pushToSynapseFxn(exprDataSTG,'MSSMRNAseq_STG.csv',MSSMSTGManifest$exprFolder@properties$id,"MSSM RNAseq STG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows")
pushToSynapseFxn(exprDataPHG,'MSSMRNAseq_PHG.csv',MSSMPHGManifest$exprFolder@properties$id,"MSSM RNAseq PHG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows")
pushToSynapseFxn(exprDataITG,'MSSMRNAseq_ITG.csv',MSSMITGManifest$exprFolder@properties$id,"MSSM RNAseq ITG expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows")
