synapseClient::synapseLogin()
summaryDf3 <- synapseClient::synTableQuery("SELECT * FROM syn10227506")@values

gs <- unique(summaryDf3$geneSet)
maxOr <- c()
medianOr <- c()
meanOr <- c()
Nam <- c()
MaxNam <- c()
MaxMethod <- c()

for (i in 1:length(unique(summaryDf3$geneSet))){
  
  In <- which(summaryDf3$geneSet == gs[i])
  In2 <- which(summaryDf3$mean_or[In]<Inf)
  if(length(In2)>0){
    maxOr <- c(maxOr,max(log(summaryDf3$mean_or[In[In2]])))
    medianOr <- c(medianOr,median(log(summaryDf3$mean_or[In[In2]])))
    meanOr <- c(meanOr,mean(log(summaryDf3$mean_or[In[In2]])))
    Nam <- c(Nam,gs[i]) 
    In3 <- which.max(log(summaryDf3$mean_or[In[In2]]))
    MaxNam <- c(MaxNam,summaryDf3$ModuleNameFull[In[In2[In3]]])
    MaxMethod <- c(MaxMethod, summaryDf3$method[In[In2[In3]]])
  }

}

DF5 <- list()
DF5$Nam <- Nam 
DF5$maxOr <- maxOr
DF5$medianOr <- medianOr
DF5$meanOr <- meanOr
DF5$MaxNam <- MaxNam
DF5$MaxMethod <- MaxMethod
DF5 <- data.frame(DF5, stringsAsFactors = F)
DF5 <- DF5[order(-DF5$medianOr),]
head(DF5, n=10)

table(DF5$MaxMethod)


