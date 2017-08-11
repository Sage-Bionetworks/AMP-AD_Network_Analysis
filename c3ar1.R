#pull bic networks;

pullAndProcess = function(gene='C3AR1'){
  #login in to synapse
  synapseClient::synapseLogin()
  #query synapse and pull down bic networks
  bicNets <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( method = 'bic' ) AND ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) )")@values
  #load bicnetworks one at a time and 1) compute degree, 2) convert to hgnc, 3) pull out c3ar
  internal = function(synId,brainRegion,gene){
    obj <- synapseClient::synGet(synId)
    load(obj@filePath)
    library(Matrix)
    foo <- as.matrix(bicNetworks$network)
    mapping <- utilityFunctions::convertEnsemblToHgnc(colnames(foo))
    foo2 <- foo+t(foo)
    degree <- rowSums(foo2)
    degreedf <- data.frame(ensembl_gene_id=names(degree),
                           degree=degree,
                           stringsAsFactors=F)
    colnames(degreedf)[2] <- paste0('degree',brainRegion)
    df <- dplyr::left_join(degreedf,
                           mapping)
    keep<-which(df$external_gene_name==gene)
    keep2 <- which(foo2[keep,]!=0)
    foo3 <- foo2[c(keep,keep2),c(keep,keep2)]
    #convert adjacency matrix to a df
    foo4 <- foo3*upper.tri(foo3)
    foo5 <- which(foo4==1,T)
    foo6 <- foo5
    foo6[,1] <- colnames(foo4)[foo5[,1]]
    foo6[,2] <- colnames(foo4)[foo5[,2]]
    colnames(foo6) <- c('ensembl_gene_id1','ensembl_gene_id2')
    foo6 <- data.frame(foo6,stringsAsFactors=F)
    foo6 <- dplyr::left_join(foo6,
                             dplyr::select(df,ensembl_gene_id,external_gene_name),
                             by=c('ensembl_gene_id1'='ensembl_gene_id'))
    colnames(foo6)[3] <- c('external_gene_name1')
    foo6 <- dplyr::left_join(foo6,
                             dplyr::select(df,ensembl_gene_id,external_gene_name),
                             by=c('ensembl_gene_id2'='ensembl_gene_id'))
    colnames(foo6)[4] <- c('external_gene_name2')
    foo6$brainRegion <- rep(brainRegion,nrow(foo6))
    foo6$external_gene_name1b <- paste0(foo6$external_gene_name1,brainRegion)
    foo6$external_gene_name2b <- paste0(foo6$external_gene_name2,brainRegion)
    return(list(degreedf=df,
                geneNghbd=foo6))
  }
  foobar <- mapply(internal,
                   bicNets$id,
                   bicNets$tissueTypeAbrv,MoreArgs = list(gene=gene),
                   SIMPLIFY=FALSE)
  #combine dfs into asingle df
  foobar2 <- lapply(foobar,function(x) x$geneNghbd)
  foobar2 <- do.call(rbind,foobar2)
  
  #
  foobar3 <- lapply(foobar,function(x) x$degreedf)
  foobar3 <- plyr::join_all(foobar3)
  model <- list()
  model$genenghbds <- foobar2
  model$degreemanifest <- foobar3
  return(model)
}
res <- pullAndProcess()
res$degreemanifest <- dplyr::mutate(res$degreemanifest,meandegree = degreeDLPFC+degreeTCX+degreeCBE+degreeFP+degreeSTG+degreePHG+degreeIFG)
degreeManifest <- res$degreemanifest
degreeManifest <- dplyr::select(degreeManifest,ensembl_gene_id,external_gene_name,degreeDLPFC,degreeTCX,degreeCBE,degreeFP,degreeSTG,degreePHG,degreeIFG)
keep1 <- c('degreeDLPFC','degreeTCX','degreeCBE','degreeFP','degreeSTG','degreePHG','degreeIFG')
degreeManifest[,keep1] <- scale(degreeManifest[,keep1])
degreeManifest <- dplyr::mutate(degreeManifest,averageDegree=degreeDLPFC+degreeTCX+degreeCBE+degreeFP+degreeSTG+degreePHG+degreeIFG)
degreeManifest$averageDegree <- degreeManifest$averageDegree/7
degreeManifest <- dplyr::arrange(degreeManifest,desc(averageDegree))
rSynapseUtilities::makeTable(degreeManifest,"metanetwork connectivity summaries",projectId = "syn2370594")

#degreeManifest1b <- scale(degreeManifest1)



write.csv(dplyr::select(res$genenghbds,external_gene_name1b,external_gene_name2b),file='~/Desktop/c3ar1.csv',quote=F)
