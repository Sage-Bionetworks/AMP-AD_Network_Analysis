synapseClient::synapseLogin()

foo <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'moduleIdentification' OR analysisType = 'consensusModuleIdentification') AND ( method = 'kmeans' OR method = 'wina' OR method = 'speakeasy' OR method = 'megena' OR method = 'wgcna' ) AND (tissueOfOrigin = 'dorsolateralPrefrontalCortex') AND (columnScaled = 'TRUE') )")

rosmapModuleManifest <- foo@values
bar <- lapply(rosmapModuleManifest$id,synapseClient::synGet)

modules <- lapply(bar,function(x){
  library(dplyr)
  synapseClient::getFileLocation(x) %>%
  data.table::fread(data.table=F)})

names(modules) <- rosmapModuleManifest$method
modules$speakEasy$V1 <- modules$kmeans$Gene.ID
colnames(modules$speakEasy) <- c('Gene.ID','Module')
colnames(modules$kmeans)[1:2] <- c('Gene.ID','Module')
modules$kmeans <- dplyr::select(modules$kmeans,-moduleLabel)

colnames(modules$megena)[1:2] <- c('Gene.ID','Module')
modules$wina <- dplyr::select(modules$wina,Geneid,module)
colnames(modules$wina)[1:2] <- c('Gene.ID','Module')

View(modules$speakEasy)
View(modules$kmeans)
View(modules$megena)
View(modules$wina)
######NMI

fxn1 <- function(x){
  y <- x$Module
  names(y) <- x$Gene.ID
  return(clue::as.cl_partition(y))
}
baz2 <- lapply(modules,fxn1)


megenaTemp <- igraph::graph_from_data_frame(modules$megena)
megenaAdj <- igraph::as_adjacency_matrix(megenaTemp)
megenaAdj <- as.matrix(megenaAdj)
moduleDefn <- unique(modules$megena$Module)
megenaAdj <- megenaAdj[-which(rownames(megenaAdj)%in%moduleDefn),moduleDefn]
notIn <- modules$kmeans$Gene.ID[which(!(modules$kmeans$Gene.ID%in%rownames(megenaAdj)))]
megenaAdj <- rbind(megenaAdj,matrix(0,length(notIn),ncol(megenaAdj)))
rownames(megenaAdj)[16716:16950] <- notIn
megenaAdj <- cbind(megenaAdj,c(rep(0,16715),rep(1,235)))
colnames(megenaAdj)[407] <- 'noMod'
megenaAdj <- megenaAdj[modules$kmeans$Gene.ID,]
foobar <- clue::as.cl_partition(as.matrix(megenaAdj))
baz2$megena <- foobar
names(baz2)[3] <- 'metanetwork'
ensembleOfCluster <- clue::cl_ensemble(list = baz2)
nmi_score <- clue::cl_agreement(ensembleOfCluster,method = "NMI")
crand_score <- clue::cl_agreement(ensembleOfCluster,method = "cRand")

listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
######overlap script
#turn into list
abcd<-lapply(modules,function(df){
  res <- lapply(unique(df$Module),listify,df$Gene.ID,df$Module)
  names(res) <- unique(df$Module)
  return(res)
})

speakeasyVmetanetworkPval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,abcd[[3]],abcd[[4]],allGenes=union(unlist(abcd[[3]]),unlist(abcd[[4]])))
speakeasyVmetanetworkOR <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,abcd[[3]],abcd[[4]],allGenes=union(unlist(abcd[[3]]),unlist(abcd[[4]])))
mv<-max(speakeasyVmetanetworkOR[-which(!is.finite(speakeasyVmetanetworkOR))])
speakeasyVmetanetworkOR[which(!is.finite(speakeasyVmetanetworkOR))] <- mv
#png(filename='~/Desktop/rosmap_msbb.png',height=4,width=6,units='in',pointsize = 5,res=300)

speakeasyVmetanetworkPval2 <- p.adjust(speakeasyVmetanetworkPval,method='fdr') %>% matrix(nrow(speakeasyVmetanetworkPval),
                                                                                          ncol(speakeasyVmetanetworkPval))


connections <- t((speakeasyVmetanetworkPval2)<0.001)
connections2 <- connections
connections2[which(connections)] <- 1
connections2[which(!connections)] <- 0
pheatmap::pheatmap(connections2*(t(speakeasyVmetanetworkOR)^(1/4)),
                   scale='none',
                   labels_col=names(abcd[[1]]),
                   xlab='ROSMAP DLPFC Modules',
                   ylab='MSBB FP Modules',
                   cluster_rows=FALSE,
                   cluster_cols=FALSE)
#dev.off()


genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList


microglialList <- GeneSets$Cell_Markers$`Zhang:Microglia`
neuralList <- GeneSets$Cell_Markers$`Zhang:Neuron`
astrocyteList <- GeneSets$Cell_Markers$`Zhang:Astrocyte`
endothelialList <- GeneSets$Cell_Markers$`Zhang:Endothelial`

library(dplyr)


adListEns <- utilityFunctions::convertHgncToEnsembl(adList)


getPvaluesAndOddsRatios <- function(list1,list2,geneSet1){
  model <- list()
  model$msbbAdPval <- lapply(list1[[1]],
                             utilityFunctions::fisherWrapperPval,
                             geneSet1,
                             unique(list2[[1]]$Gene.ID)) %>%
    unlist %>%
    p.adjust(method='fdr')
  
  model$msbbAdOR <- lapply(list1[[1]],
                           utilityFunctions::fisherWrapperOR,
                           geneSet1,
                           unique(list2[[1]]$Gene.ID)) %>%
    unlist
  
  #rosmap
  model$rosmapAdPval <- lapply(list1[[2]],
                               utilityFunctions::fisherWrapperPval,
                               geneSet1,
                               unique(list2[[2]]$Gene.ID)) %>%
    unlist %>%
    p.adjust(method='fdr')
  
  model$rosmapAdOR <- lapply(list1[[2]],
                             utilityFunctions::fisherWrapperOR,
                             geneSet1,
                             unique(list2[[2]]$Gene.ID)) %>% 
    unlist
  
  #mayo tcx
  model$mayoTcxAdPval <- lapply(list1[[3]],
                                utilityFunctions::fisherWrapperPval,
                                geneSet1,
                                unique(list2[[3]]$Gene.ID)) %>% 
    unlist %>%
    p.adjust(method='fdr')
  
  model$mayoTcxAdOR <- lapply(list1[[3]],
                              utilityFunctions::fisherWrapperOR,
                              geneSet1,
                              unique(list2[[3]]$Gene.ID)) %>% unlist
  
  return(model)
}

aaaa<-getPvaluesAndOddsRatios(abcd,modules,adListEns$ensembl_gene_id)
save(abcd,file='rosmapMods.rda')
