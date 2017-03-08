library(synapseClient)
library(githubr)
synapseLogin()



#steps: 
#1) get rid of first 4 rows
#2) split cerebellum and temporal cortex samples
#3) ensure proper location for expression data and networks on Synapse done
#4) make sure output is in proper format going forward

exprDataObj <- synGet('syn8028570')
#geneIdObj <- synGet('syn4922926')
#sampleIdObj <- synGet('syn4922923')

MayoCBEId <- 'syn7980596'
MayoTCXId <- 'syn7980592'

MayoCBEExpr <- synapseClient::Folder(name="Expression Data",parentId=MayoCBEId)
MayoCBEExpr <- synapseClient::synStore(MayoCBEExpr)
MayoCBENetworks <- synapseClient::Folder(name="Expression Networks",parentId=MayoCBEId)
MayoCBENetworks <- synapseClient::synStore(MayoCBENetworks)
metanetworkcbe <- synapseClient::Folder(name='metanetwork',parentId=MayoCBENetworks@properties$id)
megenacbe <- synapseClient::Folder(name='megena',parentId=MayoCBENetworks@properties$id)
speakeasycbe <- synapseClient::Folder(name='speakeasy',parentId=MayoCBENetworks@properties$id)
trenacbe <- synapseClient::Folder(name='trena',parentId=MayoCBENetworks@properties$id)
metanetworkcbe <- synapseClient::synStore(metanetworkcbe)
megenacbe <- synapseClient::synStore(megenacbe)
speakeasycbe <- synapseClient::synStore(speakeasycbe)
trenacbe <- synapseClient::synStore(trenacbe)
wgcnacbe <- synapseClient::synStore(wgcnacbe)


MayoTCXExpr <- synapseClient::Folder(name="Expression Data",parentId=MayoTCXId)
MayoTCXExpr <- synapseClient::synStore(MayoTCXExpr)
MayoTCXNetworks <- synapseClient::Folder(name="Expression Networks",parentId=MayoTCXId)
MayoTCXNetworks <- synapseClient::synStore(MayoTCXNetworks)
metanetworktcx <- synapseClient::Folder(name='metanetwork',parentId=MayoTCXNetworks@properties$id)
megenatcx <- synapseClient::Folder(name='megena',parentId=MayoTCXNetworks@properties$id)
speakeasytcx <- synapseClient::Folder(name='speakeasy',parentId=MayoTCXNetworks@properties$id)
trenatcx <- synapseClient::Folder(name='trena',parentId=MayoTCXNetworks@properties$id)
wgcnatcx <- synapseClient::Folder(name='wgcna',parentId=MayoTCXNetworks@properties$id)
metanetworktcx <- synapseClient::synStore(metanetworktcx)
megenatcx <- synapseClient::synStore(megenatcx)
speakeasytcx <- synapseClient::synStore(speakeasytcx)
trenatcx <- synapseClient::synStore(trenatcx)
wgcnatcx <- synapseClient::synStore(wgcnatcx)



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
tcxIds <- grep('TCX',colnames(exprData))
cerIds <- grep('CER',colnames(exprData))
exprDataTCX <- exprData[,tcxIds]
exprDataCER <- exprData[,cerIds]
#exprData <- t(exprData)
write.csv(exprDataTCX,file='MayoRNAseq_TCX.csv',quote=F)
write.csv(exprDataCER,file='MayoRNAseq_CER.csv',quote=F)

#version comment
commentTCX = "Mayo RNAseq TCX expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows"

foo = File('MayoRNAseq_TCX.csv',parentId=MayoTCXExpr@properties$id,versionComment = commentTCX)

bar2 = synGet('syn4650265',downloadFile=F)
anno = synGetAnnotations(bar2)

#annotations
synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'mayornaseq',
                               repositoryPath = 'mayoRNAseq.R')

#provenance and store
foo = synStore(foo,
               used = as.list(c('syn8028570')),
               executed = as.list(c(permLink)),
               activityName = 'Format Mayo RNAseq TCX expression Data',
               activityDescription = 'Push Mayo RNAseq TCX data into format for network pipeline')




#version comment
commentCER = "Mayo RNAseq CER expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows"

foo = File('MayoRNAseq_CER.csv',parentId=MayoCBEExpr@properties$id,versionComment = commentCER)

bar2 = synGet('syn5201007',downloadFile=F)
anno = synGetAnnotations(bar2)

#annotations
synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'mayornaseq',
                               repositoryPath = 'mayoRNAseq.R')

#provenance and store
foo = synStore(foo,
               used = as.list(c('syn8028570')),
               executed = as.list(c(permLink)),
               activityName = 'Format Mayo RNAseq CER expression Data',
               activityDescription = 'Push Mayo RNAseq CER data into format for network pipeline')

