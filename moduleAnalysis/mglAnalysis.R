
source('dataPulling/pullExpressionAndPhenoWinsorized.R')
rosmap <- geneExpressionForAnalysis$rosmapDLPFC
rosmap <- dplyr::select(rosmap,dplyr::starts_with('ENSG'))
set.seed(1)
kmeansResult <- kmeans(t(rosmap),35,nstart = 1)
df <- data.frame(cluster=as.factor(kmeansResult$cluster))
modmat <- model.matrix(~cluster-1,data = df)
foo <- MGL::MGL(data = t(rosmap),L=modmat,lambda=1)

df1 <- data.frame(gene=colnames(rosmap),module=paste0('mgl',foo$Z),stringsAsFactors=F)
df2 <- utilityFunctions::convertEnsemblToHgnc(df1$gene)
df1 <- dplyr::left_join(df1,df2,by=c('gene'='ensembl_gene_id'))
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
mglMods <- lapply(unique(df1$module),listify,df1$external_gene_name,df1$module)
names(mglMods) <- unique(df1$module)

synapseClient::synapseLogin()
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))



system.time(aaaa <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                                  mglMods,
                                                  adList2,
                                                  unique(unlist(mglMods))))
aaaa <- t(aaaa)
