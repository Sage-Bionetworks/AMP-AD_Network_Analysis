source('loadAMPADModules.R')
source('pullExpressionAndPheno.R')


allMods$ModuleNameFull <- gsub("kmeans","metanetwork",allMods$ModuleNameFull)
allMods$method <- gsub("kmeans","metanetwork",allMods$method)
allMods$ModuleName <- gsub("kmeans","metanetwork",allMods$ModuleName)

#create 7 matrices of indicator matrices
createBrainRegionMatrices <- function(x,modDefn){
  y <- dplyr::filter(modDefn,brainRegion==x)
  z <- dplyr::select(y,GeneID,method,ModuleNameFull)
  z <- tidyr::spread(z,GeneID,method)
  rownames(z) <- z$ModuleNameFull
  z <- t(z)
  #z[which(is.na(z))] <- 0
  z <- z[-1,]
  met <- unique(modDefn$method)
  z[z%in%met] <- 1
  z[is.na(z)] <- 0
  rn <- rownames(z)
  cn <- colnames(z)
  z <- matrix(as.numeric(z),nrow(z),ncol(z))
  rownames(z) <- rn
  colnames(z) <- cn
  return(z)
  #z <- data.matrix(z)
}

foo3 <- dplyr::group_by(allMods,method,Module,brainRegion) %>%dplyr::summarise(length(method))


moduleAdjacencies <- lapply(unique(allMods$brainRegion),
                            createBrainRegionMatrices,
                            allMods)
names(moduleAdjacencies) <- unique(allMods$brainRegion)

####kmeansAIC fxn
kmeansAIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

computeAICfxn <- function(k,baz2,nstartPar){
  metanetwork2 <- kmeans(baz2,k,nstart=nstartPar);
  foob<-(kmeansAIC(metanetwork2))
  cat('k:',k,'aic',foob,'\n')
  return(foob)
}

#ROSMAP
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='DLPFC')
k <- round(median(table(foo3a$method)))
#aics <- sapply(2:20,computeAICfxn,moduleAdjacencies$DLPFC,1)
#aics2 <- sapply(8:15,computeAICfxn,moduleAdjacencies$DLPFC,10)
rosmapConsensus <- kmeans(moduleAdjacencies$DLPFC,k,nstart=20)

#FP
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='FP')
k <- round(median(table(foo3a$method)))
#aics <- sapply(2:20,computeAICfxn,moduleAdjacencies$FP,1)
#aics2 <- sapply(8:15,computeAICfxn,moduleAdjacencies$FP,10)
msbbfpConsensus <- kmeans(moduleAdjacencies$FP,k,nstart=20)

#STG
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='STG')
k <- round(median(table(foo3a$method)))
msbbstgConsensus <- kmeans(moduleAdjacencies$STG,k,nstart=20)

#IFG
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='IFG')
k <- round(median(table(foo3a$method)))
msbbIFGConsensus <- kmeans(moduleAdjacencies$IFG,k,nstart=20)

#PHG
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='PHG')
k <- round(median(table(foo3a$method)))
msbbPHGConsensus <- kmeans(moduleAdjacencies$PHG,k,nstart=20)

#TCX
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='TCX')
k <- round(median(table(foo3a$method)))
msbbTCXConsensus <- kmeans(moduleAdjacencies$TCX,k,nstart=20)

#CER
set.seed(123)
foo3a<-dplyr::filter(foo3,brainRegion=='CER')
k <- round(median(table(foo3a$method)))
msbbCERConsensus <- kmeans(moduleAdjacencies$CER,k,nstart=20)

consMod <- list()
consMod$DLPFC <- rosmapConsensus$cluster
consMod$FP <- msbbfpConsensus$cluster
consMod$STG <- msbbstgConsensus$cluster
consMod$IFG <- msbbIFGConsensus$cluster
consMod$PHG <- msbbPHGConsensus$cluster
consMod$TCX <- msbbTCXConsensus$cluster
consMod$CER <- msbbCERConsensus$cluster

makeDf <- function(x){
  foo1 <- data.frame(GeneID=names(x),
                     Module=x,
                     stringsAsFactors=F)
  return(foo1)
}
consMod2 <- lapply(consMod,makeDf)
makeDfFull <- function(df,br){
  foo2 <- dplyr::mutate(df,method=rep('consensus',nrow(df)))
  foo2 <- dplyr::mutate(foo2,ModuleName=paste0(foo2$method,foo2$Module))
  geneMat <- utilityFunctions::convertEnsemblToHgnc(foo2$GeneID)
  foo2 <- dplyr::left_join(foo2,geneMat,by=c('GeneID'='ensembl_gene_id'))
  foo2 <- dplyr::mutate(foo2,brainRegion = rep(br,nrow(foo2)))
  foo2 <- dplyr::mutate(foo2,ModuleNameFull = paste0(foo2$ModuleName,foo2$brainRegion))
  return(foo2)
}
consMod3 <- mapply(makeDfFull,
                   consMod2,
                   names(consMod),
                   SIMPLIFY = F)
consMod4 <- do.call(rbind,consMod3)
allMods2 <- rbind(allMods,consMod4)
rSynapseUtilities::makeTable(allMods2,
                             tableName = "complete module manifest",
                             "syn2370594")
