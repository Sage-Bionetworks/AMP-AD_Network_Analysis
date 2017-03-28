#network comparison script
synapseClient::synapseLogin()
foo <- synapseClient::synQuery('select * from file where method==\'bic\' and projectId==\'syn2370594\'')

bar <- dplyr::select(foo,file.name,file.id,file.versionComment)

loadBic <- function(synId){
  synapseClient::synapseLogin()
  foo<-synapseClient::synGet(synId)
  load(foo@filePath)
  net <- as.matrix(bicNetworks$network)
  net <- net+t(net)
  return(net)
}

ampNetworks <- lapply(bar$file.id,loadBic)
names(ampNetworks) <- c('rosmap',
                        'mayoTcx',
                        'mayoCer',
                        'mssmFp',
                        'mssmStg',
                        'mssmPhg',
                        'mssmIfg',
                        'mcadgsTCX',
                        'mcadgsCER')
rownames(bar) <- names(ampNetworks)


nedges <- sapply(ampNetworks,function(x) sum(x))
bar$nedges <- nedges

#methodManifest <- expand.grid(names(ampNetworks),names(ampNetworks))
#methodManifest <- dplyr::filter(methodManifest,Var1!=Var2)
library(dplyr)
methodManifest <- combn(names(ampNetworks),2) %>% t
colnames(methodManifest) <- c('network1','network2')
methodManifest <- data.frame(methodManifest,
                             stringsAsFactors=F)


wrapperFxn <- function(method1,method2,allnets){
  cat('method1:',method1,'method2:',method2,'\n')
  return(metanetwork::compareTwoNetworks(allnets[[method1]],allnets[[method2]]))
}

networkCompare <- mapply(wrapperFxn,
                         methodManifest[,1],
                         methodManifest[,2],
                         MoreArgs = list(allnets=ampNetworks),
                         SIMPLIFY = FALSE)

####extract odds ratios
grabOddsRatios <- function(x){
  return(x[[1]]$estimate)
}
grabPvalue <- function(x){
  return(x[[1]]$p.value)
}
getSynIds <- function(x,networkManifest){
  return(networkManifest[x,'file.id'])
}
getNetwork1Overlap <- function(x){
  return(x[[2]])
}
getNetwork2Overlap <- function(x){
  return(x[[3]])
}
getNedges <- function(x,networkManifest){
  return(networkManifest[x,'nedges'])
}


ors <- sapply(networkCompare,grabOddsRatios)
pvals <- sapply(networkCompare,grabPvalue)
network1synId <- sapply(methodManifest$network1,getSynIds,bar)
network2synId <- sapply(methodManifest$network2,getSynIds,bar)
network1Overlap <- sapply(networkCompare,getNetwork1Overlap)
network2Overlap <- sapply(networkCompare,getNetwork2Overlap)
#colnames(methodManifest) <- c('network1','network2')
#methodManifest <- data.frame(methodManifest,
#                             stringsAsFactors = F)
methodManifest$ors <- ors
methodManifest$network1synId <- network1synId
methodManifest$network2synId <- network2synId
methodManifest$network1Overlap <- network1Overlap
methodManifest$network2Overlap <- network2Overlap
methodManifest$network1Nedges <- sapply(methodManifest$network1,getNedges,bar)
methodManifest$network2Nedges <- sapply(methodManifest$network2,getNedges,bar)

methodManifest2 <- dplyr::select(methodManifest,network1,network2,ors)
distMat <- tidyr::spread(methodManifest2,
                         network1,
                         ors)
rownames(distMat) <- distMat$network2
distMat <- dplyr::select(distMat,-network2)
missingRow <- setdiff(colnames(distMat),rownames(distMat))
missingCol <- setdiff(rownames(distMat),colnames(distMat))
distMat[,missingCol] <- NA
distMat[missingRow,] <- NA
distMat <- distMat[,rownames(distMat)]
distMat[(is.na(distMat))] <- 0
distMat2 <- distMat + t(distMat)
View(distMat2)
svd2 <- svd(distMat2)
library(wesanderson)

makePcPlot <- function(svdecomp,i,j,nameVec){
  rangex <- max(svd2$u[,i])-min(svd2$u[,i])
  rangey <- max(svd2$u[,j])-min(svd2$u[,j])
  xmin <- min(svd2$u[,i])-0.10*rangex
  xmax <- max(svd2$u[,i])+0.10*rangex
  ymin <- min(svd2$u[,j])-0.10*rangey
  ymax <- max(svd2$u[,j])+0.10*rangey
  plot(svd2$u[,i],
       svd2$u[,j],
       col='white',
       xlim=c(xmin,xmax),
       ylim=c(ymin,ymax),
       xlab = paste0('PCA',i),
       ylab = paste0('PCA',j))
  text(svd2$u[,i],
       svd2$u[,j],
       nameVec,
       col=wes_palette("Rushmore",7,type='continuous'))
}
makePcPlot(svd2,4,5,rownames(distMat))
#par(mar=c(12,4,6,6))
pheatmap::pheatmap(distMat2)


