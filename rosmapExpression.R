library(synapseClient)
library(githubr)
synapseLogin()

exprDataObj <- synGet('syn4922930')
geneIdObj <- synGet('syn4922926')
sampleIdObj <- synGet('syn4922923')

library(data.table)
exprData <- fread(exprDataObj@filePath,data.table=F)
geneId <- fread(geneIdObj@filePath,data.table=F,header = F)
View(geneId)
sampleId <- fread(sampleIdObj@filePath,data.table=F,header=F)
View(sampleId)

rownames(exprData) <- sampleId[,1]
colnames(exprData) <- geneId[,1]

write.csv(exprData,file='ROSMAP_Expression.csv',quote=F)

#version comment
comment = "ROSMAP expression data that was processed by Chris Gaiteri, with the gene names and sample ids added as row and column names"

foo = File('ROSMAP_Expression.csv',parentId='syn7981630',versionComment = comment)

bar2 = synGet('syn3505732',downloadFile=F)
anno = synGetAnnotations(bar2)

#annotations
synSetAnnotations(foo) = as.list(anno)


permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               repositoryPath = 'rosmapExpression.R')

#provenance and store
foo = synStore(foo,
               used = as.list(c('syn4922930',
                                'syn4922926',
                                'syn4922923')),
               executed = as.list())





