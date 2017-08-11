synapseClient::synapseLogin()
#mods <- list()

#mods$DLPFC <- synapseClient::synTableQuery('SELECT * FROM syn9705614')@values
#mods$IFG <- synapseClient::synTableQuery('SELECT * FROM syn9730668')@values
#mods$STG <- synapseClient::synTableQuery('SELECT * FROM syn9730669')@values
#mods$TCX <- synapseClient::synTableQuery('SELECT * FROM syn9730674')@values
#mods$CER <- synapseClient::synTableQuery('SELECT * FROM syn9730675')@values
#mods$PHG <- synapseClient::synTableQuery('SELECT * FROM syn9730672')@values
#mods$FP <- synapseClient::synTableQuery('SELECT * FROM syn9737595')@values

#get all mods
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10163855")@values

#convert to brain region based list
mods <- lapply(unique(allMods$brainRegion),
               function(x,y){
                  foo1 <- dplyr::filter(y,brainRegion==x)
                  return(foo1)},
               allMods)
names(mods) <- unique(allMods$brainRegion)


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
  baz2 <- lapply(x1,fxn1)
  
  megenaTemp <- igraph::graph_from_data_frame(x1$megena)
  megenaAdj <- igraph::as_adjacency_matrix(megenaTemp)
  megenaAdj <- as.matrix(megenaAdj)
  moduleDefn <- unique(x1$megena$Module)
  megenaAdj <- megenaAdj[-which(rownames(megenaAdj)%in%moduleDefn),moduleDefn]
  notIn <- x1$metanetwork$GeneID[which(!(x1$metanetwork$GeneID%in%rownames(megenaAdj)))]
  ind3 <- nrow(megenaAdj)
  megenaAdj <- rbind(megenaAdj,matrix(0,length(notIn),ncol(megenaAdj)))
  ind1 <- nrow(x1$metanetwork) - length(notIn)+1
  ind2 <- nrow(x1$metanetwork)
  
  ind4 <- length(notIn)
  rownames(megenaAdj)[ind1:ind2] <- notIn
  megenaAdj <- cbind(megenaAdj,c(rep(0,ind3),rep(1,ind4)))
  colnames(megenaAdj)[ncol(megenaAdj)] <- 'noMod'
  megenaAdj <- megenaAdj[x1$metanetwork$GeneID,]
  foobar <- clue::as.cl_partition(as.matrix(megenaAdj))
  baz2$megena <- foobar
  #names(baz2)[3] <- 'metanetwork'
  ensembleOfCluster <- clue::cl_ensemble(list = baz2)
  return(ensembleOfCluster)
}


nmi <- lapply(mods[1],getClueMods)
nmi_score <- lapply(nmi,clue::cl_agreement,method='NMI')

convertNmiToFlatTable <- function(x){
  x <- as.matrix(x)
  #print(x)
  foo1 <- which(lower.tri(x),T)
  foo2 <- c(as.matrix(x))[which(lower.tri(x))]
  foo3 <- foo1
  foo3[,1] <- rownames(x)[foo1[,1]]
  foo3[,2] <- colnames(x)[foo1[,2]]
  df1 <- cbind(foo3,foo2)
  colnames(df1) <- c('method1','method2','nmi')
  df1 <- data.frame(df1,stringsAsFactors=F)
  return(df1)
}
nmi_score_lf_list <- lapply(nmi_score,convertNmiToFlatTable)

addBrainRegion <- function(x,y){
  x <-cbind(x,rep(y,nrow(x)))
  colnames(x)[ncol(x)] <- 'tissue'
  return(x)
}
nmi_score_lf_list <- mapply(addBrainRegion,nmi_score_lf_list,names(nmi_score_lf_list),SIMPLIFY=F)
nmi_score_lf_list <- do.call(rbind,nmi_score_lf_list)
View(nmi_score_lf_list)
rSynapseUtilities::makeTable(nmi_score_lf_list,'pairwise nmi scores with consensus','syn2370594')

png(file='~/Desktop/nmirosmap.png',
    height=800,
    width=1000,
    res=120,
    pointsize = 30)
pheatmap::pheatmap(data.matrix(nmi_score[[1]]),fontsize = 25)
dev.off()
