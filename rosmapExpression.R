library(synapseClient)
library(githubr)
synapseLogin()

exprDataObj <- synGet('syn8018356')
#geneIdObj <- synGet('syn4922926')
#sampleIdObj <- synGet('syn4922923')

library(data.table)
exprData <- fread(exprDataObj@filePath,data.table=F)
#geneId <- fread(geneIdObj@filePath,data.table=F,header = F)
#View(geneId)
#sampleId <- fread(sampleIdObj@filePath,data.table=F,header=F)
#View(sampleId)

#rownames(exprData) <- sampleId[,1]
#colnames(exprData) <- geneId[,1]
geneKey <- exprData[,1:3]
#rownames(exprData) <- exprData[,1]
exprData <- exprData[,-c(1:3)]
exprData <- t(exprData)
colnames(exprData) <- geneKey$ensembl_gene_id
exprData <- t(exprData)
write.csv(exprData,file='ROSMAP_Expression.csv',quote=F)

#version comment
comment = "ROSMAP expression data processed with the RNAseq reprocessing pipeline in the AMP-AD consortia, residualized for batch and technical confounds with genes as rows"

foo = File('ROSMAP_Expression.csv',parentId='syn7981630',versionComment = comment)

bar2 = synGet('syn3505732',downloadFile=F)
anno = synGetAnnotations(bar2)

#annotations
synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'rosmap',
                               repositoryPath = 'rosmapExpression.R')

#provenance and store
foo = synStore(foo,
               used = as.list(c('syn8018356')),
               executed = as.list(c(permLink)),
               activityName = 'Format ROSMAP expression Data',
               activityDescription = 'Push ROSMAP data into format for network pipeline')





