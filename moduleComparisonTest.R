synapseClient::synapseLogin()
foo <- synapseClient::synQuery('select name,id,method from File where projectId==\'syn2370594\' and analysisType==\'moduleIdentification\' and study==\'ROSMAP\'')
bar <- lapply(foo$File.id,synapseClient::synGet)
names(bar) <- foo$File.method
library(dplyr)
baz <- lapply(bar,synapseClient::getFileLocation) %>%
  lapply(data.table::fread,data.table=F)

fxn1 <- function(x){
  y <- x$moduleLabel
  names(y) <- x$Gene.ID
  return(clue::as.cl_partition(y))
}
baz2 <- lapply(baz,fxn1)

ensembleOfCluster <- clue::cl_ensemble(list = baz2)
cl_agreement(ensembleOfCluster,method = "euclidean")
bar22 <- cl_consensus(ensembleOfCluster)
