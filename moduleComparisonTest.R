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
clue::cl_agreement(ensembleOfCluster,method = "NMI")
bar22 <- clue::cl_consensus(ensembleOfCluster)
moduleAssignments <- apply(bar22$.Data,1,function(x){if(sum(x>.5)>0){return(which(x>.5))}else{return(NA)}})
v1 <-utilityFunctions::convertEnsemblToHgnc(names(moduleAssignments))
moduleAssignments2 <- data.frame(moduleLabels=moduleAssignments,
                                 ensembl_gene_id=names(moduleAssignments),stringsAsFactors=F)
v1 <- dplyr::left_join(v1,moduleAssignments2,"ensembl_gene_id")
v_29 <- dplyr::filter(v1,moduleLabels==169)
cat(v_29$external_gene_name,file='~/Desktop/mod29.csv',sep='\n')


#########overlaps

