load('aggregate_module_mainfest.rda')

#CBE Turquoise, PHG Turquoise, STG Turquoise, TCX turquoise, IFG Turquoise, DLPFC blue, TCX blue

mods <- c("aggregateCBEturquoiseCBE",
          "aggregatePHGturquoisePHG",
          "aggregateSTGturquoiseSTG",
          "aggregateTCXturquoiseTCX",
          "aggregateIFGturquoiseIFG",
          "aggregateDLPFCblueDLPFC",
          "aggregateTCXblueTCX")
manifest <- fullManifest[mods]

edgeLists <- lapply(manifest,function(x){
  library(Matrix)
  foobar <- utilityFunctions::convertAdjacencyToEdgeList(x$adjacencyMatrix)
  foobar1 <- utilityFunctions::convertEnsemblToHgnc(unique(c(foobar)))
  foobar <- data.frame(foobar,stringsAsFactors=F)
  colnames(foobar) <- c('node1','node2')
  foobar <- dplyr::left_join(foobar,foobar1,by=c('node1'='ensembl_gene_id'))
  foobar <- dplyr::left_join(foobar,foobar1,by=c('node2'='ensembl_gene_id'))
  return(foobar)
})

filenames <- c('cbe_turquoise.csv',
               'phg_turquoise.csv',
               'stg_turquoise.csv',
               'tcx_turquoise.csv',
               'ifg_turquoise.csv',
               'dlpfc_turquoise.csv',
               'tcx_blue.csv')
mapply(function(x,y){
  write.csv(x[,3:4],file=y,quote=F)
},edgeLists,filenames)

#make table of cell types
synapseClient::synapseLogin()
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
cellTypeDf <- utilityFunctions::list2df(GeneSets$Cell_Markers)
write.csv(cellTypeDf,file='cellTypeAnnotations.csv',quote=F)
#GeneSets$Cell_Markers
