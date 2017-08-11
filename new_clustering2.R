synapseClient::synapseLogin()
modDefs <- synapseClient::synTableQuery("select * from syn10146524")@values
head(modDefs)


modDefsDLPFC <- dplyr::filter(modDefs,brainRegion == 'DLPFC')
convertFormat <- function(methodVal,modTab){
  foo <- dplyr::filter(modTab,method==methodVal)
  bar <- dplyr::select(foo,GeneID,Module)
  colnames(bar) <- c('Gene.ID','moduleNumber')
  return(bar)
}

partition.adj <- lapply(unique(modDefsDLPFC$method),convertFormat,modDefsDLPFC)
names(partition.adj) <- unique(modDefsDLPFC$method)
library(dplyr)
partition.adj <- mapply(function(mod, method){
  mod = mod %>%
    dplyr::select(Gene.ID, moduleNumber) %>%
    dplyr::mutate(value = 1,
                  moduleNumber = paste0(method,'.',moduleNumber)) %>%
    tidyr::spread(moduleNumber, value)},
  partition.adj,
  names(partition.adj),
  SIMPLIFY = F) %>%
  plyr::join_all(type="full")
partition.adj[is.na(partition.adj)] <- 0
rownames(partition.adj) <- partition.adj$Gene.ID
partition.adj$Gene.ID <- NULL
set.seed(1)
# Randomise gene order
partition.adj <- partition.adj[sample(1:dim(partition.adj)[1], dim(partition.adj)[1]), ]
dim(partition.adj)

# pull pairwise
pairwiseObj <- synapseClient::synGet("syn10146717")
foo2 <- data.table::fread(pairwiseObj@filePath,data.table=F)
dlpfc <- intersect(grep('DLPFC',foo2$ModuleNameFull),grep('DLPFC',foo2$category))
foo3 <- foo2[dlpfc,]
library(dplyr)
foo4 <- dplyr::select(foo3,ModuleNameFull,category,fisherPval) %>%
  tidyr::spread(category,fisherPval)
rownames(foo4) <- foo4$ModuleNameFull
foo4 <- dplyr::select(foo4,-ModuleNameFull)
foo4 <- data.matrix(foo4)
foo5 <- foo4
nsig<-sum(p.adjust(c(foo4[which(upper.tri(foo4))]),method='fdr')< 0.05)
pval<-sort(c(foo4[which(upper.tri(foo4))]))[nsig]

foo5[which(foo4 < pval)] <- 1
foo5[which(foo4 >=pval)] <- 0

foo6 <- dplyr::select(foo3,ModuleNameFull,category,fisherOR) %>%
  tidyr::spread(category,fisherOR)
rownames(foo6) <- foo6$ModuleNameFull
foo6 <- dplyr::select(foo6,-ModuleNameFull)
foo6 <- data.matrix(foo6)

foo7 <- foo6
foo7[which(!is.finite(foo6))]<- NA
mv <- max(foo7,na.rm=T)
foo7[which(is.na(foo7))] <- mv

foo7 <- foo7*foo5

foo7 <- foo7 + 1e-16
foo8 <-pheatmap::pheatmap(log(foo7),scale='none')



foo9 <-cutree(foo8$tree_col,k = 200)
moduleCluster <- data.frame(ModuleNameFull = names(foo9),metaCluster = foo9,stringsAsFactors=F)
moduleCluster <- dplyr::arrange(moduleCluster,metaCluster)
buildMetaClusters <- function(i,moduleCluster,modDefsDLPFC){
  ab<-dplyr::filter(moduleCluster,metaCluster==i)
  cd <- dplyr::left_join(ab,modDefsDLPFC)
  #cd <- cd[which(!duplicated(cd)),]
  return(cd)
}
afc <- lapply(unique(moduleCluster$metaCluster),
              buildMetaClusters,
              moduleCluster,
              modDefsDLPFC)
modsize <- sapply(afc,function(x){return(length(unique(x$GeneID)))})

geneLists <- lapply(afc,function(x) unique(x$external_gene_name))

genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))




system.time(aaaa <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                                  geneLists,
                                                  GeneSets$KEGG_2015,
                                                  unique(unlist(geneLists))))
