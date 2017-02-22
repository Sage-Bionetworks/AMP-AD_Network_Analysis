library(synapseClient)
synapseLogin()
winsorizeAndScale <- function(x){
  library(dplyr)
  library(utilityFunctions)
  x <- t(x)
  x <- apply(x,2,utilityFunctions::winsorize)
  x <- scale(x)
  return(x)
}
#######query private project and pull data
#dataType=mRNA
foo <- synQuery('select name,id from file where projectId==\'syn2370594\' and dataType==\'mRNA\'')

#get data
bar <- lapply(foo$file.id,synGet)

#function to load data quickly
loadData <- function(x){
  library(data.table)
  testData <- data.table::fread(x@filePath,data.table=F)
  rownames(testData) <- testData[,1]
  testData <- testData[,-1]
  return(testData)
}

#load all 7 data sets into a list
listOfData <- lapply(bar,loadData)

#winsorize and scale
reformattedData <- lapply(listOfData,winsorizeAndScale)

library(dplyr)

fileNames <- paste0('Scaled_Winsorized_',foo$file.name)
parentIds <- sapply(bar,function(x){return(x@properties$parentId)})
annotationss <- sapply(bar,function(x){return(synGetAnnotations(x))})
comments <- rep('Winsorizing and scaling',7)
useds <- sapply(bar,function(x){return(x@properties$id)})

library(githubr)

permLink1 =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'mayornaseq',
                               repositoryPath = 'makeDataForNetworkWG.R')

permLink2 =githubr::getPermlink(repository = 'Sage-Bionetworks/rSynapseUtilities',
                                ref = 'branch',
                                refName = 'dev',
                                repositoryPath = 'pushToSynapseWrapper.R')

permLink3 =githubr::getPermlink(repository = 'blogsdon/utilityFunctions',
                                ref = 'branch',
                                refName = 'master',
                                repositoryPath = 'winsorize.R')



executeds <- c(permLink1,permLink2,permLink3) %>% rep(7)

activityNames <- rep('Post Process AMP-AD RNAseq data',7)

activityDescriptions <- rep('Transpose, Winsorize, then Scale each data-set',7)

#push back to synapse
synObjs<-mapply(rSynapseUtilities::pushToSynapseWrapper,
       df=reformattedData,
       fileName=fileNames,
       synapseFolderId=parentIds,
       annos=annotations,
       comment=comments,
       usedVector=useds,
       executedVector=executeds,
       activityName1=activityNames,
       activityDescription1=activityDescriptions)

