#network comparison script
synapseClient::synapseLogin()
foo <- synapseClient::synQuery('select * from file where method==\'bic\' and projectId==\'syn2370594\'')

bar <- dplyr::select(foo,file.name,file.id,file.versionComment)

loadBic <- function(synId){
  synapseClient::synapseLogin()
  foo<-synapseClient::synGet(synId)
  load(foo@filePath)
  return(bicNetworks$network)
}

ampNetworks <- lapply(bar$file.id,loadBic)
names(ampNetworks) <- c('rosmap','mayoTcx','mayoCer','mssmFp','mssmStg','mssmPhg','mssmIfg')

#methodManifest <- expand.grid(names(ampNetworks),names(ampNetworks))
#methodManifest <- dplyr::filter(methodManifest,Var1!=Var2)
library(dplyr)
methodManifest <- combn(names(ampNetworks),2) %>% t

wrapperFxn <- function(method1,method2,allnets){
  return(metanetwork::compareTwoNetworks(allnets[[method1]],allnets[[method2]]))
}

networkCompare <- mapply(wrapperFxn,
                         methodManifest[,1],
                         methodManifest[,2],
                         MoreArgs = list(allnets=ampNetworks),
                         SIMPLIFY = FALSE)
