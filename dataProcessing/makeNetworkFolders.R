require(synapseClient)
synapseLogin()
foo <- synQuery('select name,id,parentId from folder where name==\'metanetwork\' and projectId==\'syn2370594\'')

createFolder <- function(synId){
  bar <- synapseClient::Folder(name='networks',
                               parentId=synId)
  bar <- synapseClient::synStore(bar)
  return(bar)
}

moveFile <- function(fileId,parentId){
  myFile <- synapseClient::synGet(fileId,downloadFile = FALSE);
  myFile$properties$parentId <- parentId;
  myFile <- synapseClient::synStore(myFile,forceVersion=F);
}


foo1 <- sapply(foo$folder.id,createFolder)
foo2 <- sapply(foo1,function(x){return(x@properties$id)})
foo$newFolder <-foo2

networkManifest <- synQuery('select name,id,parentId from file where analysisType==\'statisticalNetworkReconstruction\' and projectId==\'syn2370594\'')

networkManifest2 <- dplyr::left_join(networkManifest,foo,by=c("file.parentId"="folder.id"))
bar4<-mapply(moveFile,networkManifest2$file.id,networkManifest2$newFolder)

