cat('logging into Synapse...\n')
synapseClient::synapseLogin()
#grab module definitions
cat('pulling modules...\n')
allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ","syn10163855"))@values



mods <- lapply(unique(allMods$brainRegion),
               function(x,y){
                 foo1 <- dplyr::filter(y,brainRegion==x)
                 return(foo1)},
               allMods)
names(mods) <- unique(allMods$brainRegion)

keepGene <- unique(mods[[1]]$GeneID)
for (i in 2:7){
  keepGene <- intersect(keepGene,mods[[i]]$GeneID)
}
allMods2 <- dplyr::filter(allMods,GeneID%in%keepGene)
mods <- lapply(unique(allMods2$brainRegion),
               function(x,y){
                 foo1 <- dplyr::filter(y,brainRegion==x)
                 return(foo1)},
               allMods2)
names(mods) <- unique(allMods2$brainRegion)

####compute NMI

######NMI
getClueMods <- function(x){
  fxn1 <- function(x){
    y <- x$Module
    names(y) <- x$GeneID
    return(clue::as.cl_partition(y))
  }
  listify <- function(y,x){
    return(dplyr::filter(x,method==y))
  }
  x1 <- lapply(unique(x$method),listify,x)
  names(x1) <- unique(x$method)
  aa<-which(names(x1)=='metanetwork')
  if(length(aa)>0){
    names(x1)[aa] <- 'kmeans'
  }
  baz2 <- lapply(x1,fxn1)
  
  megenaTemp <- igraph::graph_from_data_frame(x1$megena)
  megenaAdj <- igraph::as_adjacency_matrix(megenaTemp)
  megenaAdj <- as.matrix(megenaAdj)
  moduleDefn <- unique(x1$megena$Module)
  megenaAdj <- megenaAdj[-which(rownames(megenaAdj)%in%moduleDefn),moduleDefn]
  notIn <- x1$kmeans$GeneID[which(!(x1$kmeans$GeneID%in%rownames(megenaAdj)))]
  ind3 <- nrow(megenaAdj)
  megenaAdj <- rbind(megenaAdj,matrix(0,length(notIn),ncol(megenaAdj)))
  ind1 <- nrow(x1$kmeans) - length(notIn)+1
  ind2 <- nrow(x1$kmeans)
  
  ind4 <- length(notIn)
  rownames(megenaAdj)[ind1:ind2] <- notIn
  megenaAdj <- cbind(megenaAdj,c(rep(0,ind3),rep(1,ind4)))
  colnames(megenaAdj)[ncol(megenaAdj)] <- 'noMod'
  megenaAdj <- megenaAdj[x1$kmeans$GeneID,]
  foobar <- clue::as.cl_partition(as.matrix(megenaAdj))
  baz2$megena <- foobar
  #names(baz2)[3] <- 'kmeans'
  return(baz2)
}

clustLists <- lapply(mods,getClueMods)
nmi_speakeasy <- lapply(clustLists,function(x) x$consensus)
nmi_speakeasy <- clue::cl_ensemble(list = nmi_speakeasy)
nmi_score <- clue::cl_agreement(nmi_speakeasy,method='NMI')

pheatmap::pheatmap(data.matrix(nmi_score),main = 'consensus')
