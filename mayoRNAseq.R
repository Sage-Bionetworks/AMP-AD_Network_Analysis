synapseClient::synapseLogin()



#steps: 
#1) get rid of first 4 rows
#2) split cerebellum and temporal cortex samples
#3) ensure proper location for expression data and networks on Synapse done
#4) make sure output is in proper format going forward

exprDataObj <- synapseClient::synGet('syn8466826')
#geneIdObj <- synGet('syn4922926')
#sampleIdObj <- synGet('syn4922923')


exprData <- data.table::fread(exprDataObj@filePath,data.table=F)
geneKey <- exprData[,1]
exprData <- exprData[,-c(1)]
exprData <- t(exprData)
colnames(exprData) <- geneKey
exprData <- t(exprData)
tcxIds <- grep('TCX',colnames(exprData))
cerIds <- grep('CER',colnames(exprData))
exprDataTCX <- exprData[,tcxIds]
exprDataCER <- exprData[,cerIds]
#exprData <- t(exprData)
write.csv(exprDataTCX,file='MayoRNAseq_TCX.csv',quote=F)
write.csv(exprDataCER,file='MayoRNAseq_CER.csv',quote=F)

#version comment
commentTCX = "Mayo RNAseq TCX expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia for network analysis"

foo = synapseClient::File('MayoRNAseq_TCX.csv',parentId='syn8257429',versionComment = commentTCX)

bar2 = synapseClient::synGet('syn4650265',downloadFile=F)
anno = synapseClient::synGetAnnotations(bar2)

#annotations
synapseClient::synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'mayo-patch-1',
                               repositoryPath = 'mayoRNAseq.R')

#provenance and store
foo = synapseClient::synStore(foo,
               used = as.list(c('syn8466826')),
               executed = as.list(c(permLink)),
               activityName = 'Format Mayo RNAseq TCX expression Data',
               activityDescription = 'Push Mayo RNAseq TCX data into format for network pipeline')




#version comment
commentCER = "Mayo RNAseq CER expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes in the correct format for metanetwork as rows"

foo = synapseClient::File('MayoRNAseq_CER.csv',parentId='syn8257427',versionComment = commentCER)

bar2 = synapseClient::synGet('syn5201007',downloadFile=F)
anno = synapseClient::synGetAnnotations(bar2)

#annotations
synapseClient::synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'mayo-patch-1',
                               repositoryPath = 'mayoRNAseq.R')

#provenance and store
foo = synapseClient::synStore(foo,
               used = as.list(c('syn8466826')),
               executed = as.list(c(permLink)),
               activityName = 'Format Mayo RNAseq CER expression Data',
               activityDescription = 'Push Mayo RNAseq CER data into format for network pipeline')

