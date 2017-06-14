library(synapseClient)
library(githubr)
synapseLogin()



#steps: 
#1) get rid of first 4 rows
#2) split cerebellum and temporal cortex samples
#3) ensure proper location for expression data and networks on Synapse done
#4) make sure output is in proper format going forward

cerExprDataObj <- synGet('syn5605698')
tcxExprDataObj <- synGet('syn5605688')
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


MCADGS_CER_Id <- 'syn7980585'
MCADGS_TCX_Id <- 'syn7980579'

MCADGS_CER_Manifest <- populateFolders(MCADGS_CER_Id)
MCADGS_TCX_Manifest <- populateFolders(MCADGS_TCX_Id)


library(data.table)

####temporal cortex

collapseExprToGene <- function(exprData){
  geneKey <- exprData[,1:2]
  exprData <- exprData[,-c(1:2)]
  exprData <- t(exprData)
  colnames(exprData) <- geneKey$illumina_humanht_12_v4
  
  library(biomaRt)
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='www.ensembl.org')
  
  geneKey2 <- utilityFunctions::convertIlluminaHumanHt12v4ToEnsembl(geneKey$illumina_humanht_12_v4)
  
  library(dplyr)
  geneKey3 <- dplyr::left_join(geneKey,geneKey2,by="illumina_humanht_12_v4")
  
  ####remove all illumina ids that map to multiple ensembl ids
  getEnsemblIds <- function(illuminaId,geneKey3){
    foo1 <- dplyr::filter(geneKey3,illumina_humanht_12_v4==illuminaId)
    return(unique(foo1$ensembl_gene_id))
  }
  illuminaEnsMap <- lapply(unique(geneKey3$illumina_humanht_12_v4),getEnsemblIds,geneKey3)
  names(illuminaEnsMap) <- unique(geneKey3$illumina_humanht_12_v4)
  lenVec <- sapply(illuminaEnsMap,length)
  keepIlluminaIds <- names(illuminaEnsMap)[which(lenVec==1)]
  geneKey4 <- dplyr::filter(geneKey3,illumina_humanht_12_v4%in%keepIlluminaIds)
  
  exprData <- exprData[,which(colnames(exprData)%in%keepIlluminaIds)]
  
  #####get mapping between unique ensembl ids and illumina ids
  getIlluminaIds <- function(ensemblId,geneKey4){
    foo1 <- dplyr::filter(geneKey4,ensembl_gene_id==ensemblId)
    return(unique(foo1$illumina_humanht_12_v4))
  }
  ensemblIllMap <- lapply(unique(geneKey4$ensembl_gene_id),getIlluminaIds,geneKey4)
  names(ensemblIllMap) <- unique(geneKey4$ensembl_gene_id)
  lenVec2 <- sapply(ensemblIllMap,length)
  ensemblIllMap <- ensemblIllMap[-which(lenVec2==0)]
  lenVec2 <- sapply(ensemblIllMap,length)
  
  
  #average expression across probes that map to same gene
  produceGeneLevelMatrix <- function(illId,expressionData){
    library(dplyr)
    res <- expressionData[,illId] %>%
      as.matrix %>%
      rowMeans
    return(res)
  }
  
  #make into a data frame and write to file
  listVersion <- lapply(ensemblIllMap,
                        produceGeneLevelMatrix,
                        exprData)
  
  
  names(listVersion) <- names(ensemblIllMap)
  tcxExprData <- as.data.frame(listVersion,stringsAsFactors=F)
  tcxExprData <- t(tcxExprData)
  #write.csv(tcxExprData,file=fileName,quote=F)
  return(tcxExprData)
}
exprData <- fread(tcxExprDataObj@filePath,data.table=F)
exprDataTCX<-collapseExprToGene(exprData)
exprData <- fread(cerExprDataObj@filePath,data.table=F)
exprDataCER<-collapseExprToGene(exprData)

pushToSynapseFxn <- function(exprMat,fileName,parentId,comment,annoId,usedId){
  write.csv(exprMat,file=fileName,quote=F)
  foo = File(fileName,parentId=parentId,versionComment = comment)
  bar2 = synGet(annoId,downloadFile=F)
  anno = synGetAnnotations(bar2)
  synSetAnnotations(foo) = as.list(anno)
  permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                                 ref = 'branch',
                                 refName = 'mcadgs',
                                 repositoryPath = 'mcadgsExpr.R')
  foo = synStore(foo,
                 used = as.list(c(usedId)),
                 executed = as.list(c(permLink)),
                 activityName = 'Format MCADGS Array Expression Data',
                 activityDescription = 'Push MCADGS Array Expression Data into format for network pipeline')
}

pushToSynapseFxn(exprDataCER,"MCADGS_CER.csv",MCADGS_CER_Manifest$exprFolder@properties$id,"MCADGS Cerebellum array expression data processed in the AMP-AD consortia with probe ids that map to multiple genes discarded and probes that mapped to the same gene averaged","syn3256501","syn5605698")
pushToSynapseFxn(exprDataTCX,"MCADGS_TCX.csv",MCADGS_TCX_Manifest$exprFolder@properties$id,"MCADGS temporal cortex expression data processed in the AMP-AD consortia with probe ids that map to multiple genes discarded and probes that mapped to the same gene averaged","syn3617054","syn5605688")

