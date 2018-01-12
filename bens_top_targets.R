synapser::synLogin()

#pull modules from synapse
aggregateModules <- synapser::synTableQuery("select * from syn11182793")
aggMods <-aggregateModules$asDataFrame()

#pull network from synapse
foo <- synapser::synTableQuery('select * from syn8681664 where method=\'bic\' and assay =\'RNAseq\'')
foo <- foo$asDataFrame()

loadBic <- function(synId){
  foo<-synapser::synGet(synId)
  load(foo$path)
  net <- as.matrix(bicNetworks$network)
  net <- net+t(net)
  return(net)
}



ampNetworks <- lapply(foo$id,loadBic)
names(ampNetworks) <- foo$tissueTypeAbrv

rownames(foo) <- names(ampNetworks)


nedges <- sapply(ampNetworks,function(x) sum(x))
foo$nedges <- nedges

rm(ampNetworks)
gc()
networkGraphs <- lapply(ampNetworks,
                        igraph::graph_from_adjacency_matrix,
                        mode='undirected')
#rosmap <- igraph::graph_from_adjacency_matrix(bicNetworks$network,
#                                              mode='undirected')


collateGraphStatistics <- function(graph){
  model <- list()
  cat('Computing Degree...\n')
  model$degree <- igraph::degree(graph)
  #model$alpha_centrality <- igraph::alpha_centrality(graph,tol=1e-14)
  cat('Computing Authority Score...\n')
  model$authority_score <- igraph::authority_score(graph)$vector
  
  cat('Computing Closeness...\n')
  model$closeness <- igraph::closeness(graph)
  
  cat('Computing Eccentricity...\n')
  model$eccentricity <- igraph::eccentricity(graph)
  
  cat('Computing Eigenvector Centrality...\n')
  model$eigen_centrality <- igraph::eigen_centrality(graph)$vector
  
  cat('Computing Betweeness Centrality...\n')
  model$centr_betw <- igraph::betweenness(graph)
  
  cat('Computing page rank...\n')
  model$pagerank <- igraph::page.rank(graph)$vector
  
  cat('Computing Transitivity...\n')
  model$transitivity <- igraph::transitivity(graph,
                                             type='undirected')
  
  
  return(model)
}

networkProperties <- lapply(networkGraphs,collateGraphStatistics)
str(networkProperties)

networkProperties2 <- lapply(networkProperties,
                             function(x){y <- data.frame(x[-which(names(x)=='transitivity')],stringsAsFactors=F);y$GeneID <- names(x$degree);return(y);})

dummyFun <- function(x,y){
  merge(x,y,by='GeneID')
}
                                          
networkProperties3 <- mapply(function(x,y){
  colnames(y)[1:7] <- paste0(x,colnames(y)[1:7]);return(y)},
  names(networkProperties2),
  networkProperties2,SIMPLIFY = F)


networkProperties4 <- Reduce(merge,networkProperties3)
networkProperties5 <- dplyr::select(networkProperties4,-GeneID)


genesets1 <- synapser::synGet('syn5923958')
load(genesets1$path)
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]

adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))
adList2$ad_gwas <- c(adList2$ad_gwas,'ABI3','PLCG2')
adList2$ad_gwas <- c(adList2$ad_gwas,'SPI1')

fb<-readLines('ng_gene')
adList2$ad_gwas <- union(adList2$ad_gwas,fb)

convertAd <- utilityFunctions::convertHgncToEnsembl(adList2$ad_gwas)
adGene <- networkProperties4$GeneID%in% convertAd$ensembl_gene_id



#dv <- rowMeans(dplyr::select(networkProperties5,dplyr::ends_with('degree')))
#ev <- rowMeans(dplyr::select(networkProperties5,dplyr::ends_with('authority_score')))
#tcx <- dplyr::select(networkProperties5,dplyr::starts_with('TCX'))
#tcx <- dplyr::select(tcx,-TCXeigen_centrality)
#summary(glm(adGene ~ .,data=tcx))
fun <- function(x,y){
  return(summary(glm(y~x,family='binomial'))$coef[2,4])
}
fun2 <- function(x,y){
  return(summary(glm(y~x,family='binomial'))$coef[2,1])
}
networkProperties6 <- dplyr::select(networkProperties5,-dplyr::ends_with('eigen_centrality'))
networkProperties6 <- scale(networkProperties6)


####build deg set
degResObj <- synapser::synGet("syn10496554")
load(degResObj$path)

foo <- utilityFunctions::list2df(amp.ad.de.geneSets)
View(foo)

#get rid of missing values
foo <- dplyr::filter(foo,value!='')
View(foo)

#Add ENSEMBL Ids back
map <- utilityFunctions::convertHgncToEnsembl(foo$value)
foo <- dplyr::left_join(foo,map,by=c('value'='external_gene_name'))
foo <- foo[!duplicated(foo),]
foo2 <- tidyr::drop_na(foo)
foo2 <- dplyr::select(foo2,key,ensembl_gene_id)
foo2$value <- rep(1,nrow(foo2))
foo3 <- tidyr::spread(foo2,key,value,fill=0)

combinedFeatureSet <- dplyr::left_join(networkProperties4,foo3,by=c('GeneID'='ensembl_gene_id'))
combinedFeatureSet2 <- dplyr::select(combinedFeatureSet,-dplyr::ends_with('eigen_centrality'))
combinedFeatureSet2 <- dplyr::select(combinedFeatureSet2,-GeneID)
combinedFeatureSet2[is.na(combinedFeatureSet2)] <- 0
adGene2 <- combinedFeatureSet$GeneID %in% convertAd$ensembl_gene_id
wZ <- which(colSums(combinedFeatureSet2)<3 & (1:ncol(combinedFeatureSet2))%in%(43:107))
combinedFeatureSet2 <- combinedFeatureSet2[,-wZ]
combinedFeatureSet2 <- scale(combinedFeatureSet2)

res_cf <- apply(combinedFeatureSet2,2,fun,adGene2)
lasso <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial')
#####run analyses

res<-apply(networkProperties6,2,fun,adGene)
res2<-apply(networkProperties6,2,fun2,adGene)

gap::qqunif(res,main='network features AD association')

lasso <- glmnet::cv.glmnet(y=adGene,x=networkProperties6,family='binomial')
lasso2 <- glmnet::glmnet(y=adGene,x=networkProperties6,family='binomial')
betaScore <- lasso$glmnet.fit$beta[,which(lasso$lambda==lasso$lambda.min)]
score3 <- networkProperties6[,names(betaScore)]%*%betaScore
names(score3) <- networkProperties4$GeneID
score3Df <- data.frame(gene=names(score3),score5=score3,stringsAsFactors = F)
mapTable<-dplyr::left_join(mapTable,score3Df,by=c('ensembl_gene_id'='gene'))
mapTable$positiveControl <- mapTable$external_gene_name %in% adList2$ad_gwas
mapTable<-dplyr::arrange(mapTable,desc(score4))
write.csv(dplyr::select(mapTable,external_gene_name,positiveControl,score4),file='driverScore.csv',quote=F,row.names=F)

View(mapTable)


res_sig <- p.adjust(res,method='bonferroni')
res_sig <- res_sig[res_sig<=0.05]

res2_sig <- res2[names(res_sig)]
networkProperties7 <- networkProperties6[,names(res_sig)]
masterRank <- networkProperties7%*%res2_sig
names(masterRank) <- networkProperties4$GeneID
sort(masterRank,decreasing=T)[1:5]
mapTable <- utilityFunctions::convertEnsemblToHgnc(names(masterRank))

masterRankDf <- data.frame(gene = names(masterRank),score=masterRank,stringsAsFactors = F)
mapTable <- dplyr::left_join(mapTable,masterRankDf,by=c('ensembl_gene_id'='gene'))

networkProperties8 <- scale(networkProperties5)
null = glm(adGene~1,data=networkProperties5,family='binomial')

full=glm(adGene~., data=networkProperties5,family='binomial')
res<-step(null, scope=list(lower=null, upper=full), direction="forward",steps=8)

score2 <- networkProperties5[,rownames(summary(foo)$coef[2:7,])]%*%summary(foo)$coef[2:7,1]
names(score2) <- networkProperties4$GeneID
score2Merge <- data.frame(gene=names(score2),score2=score2,stringsAsFactors=F)
mapTable <- dplyr::left_join(mapTable,score2Merge,by=c('ensembl_gene_id'='gene'))

View(mapTable)

####to do before call tomorrow
#1. make a slide describing what I did (networks -> features -> ad gene association score)
#2. make a slide on the qq plots
#3. make a slide on the positie controls (SPI1, ABI3, other genes etc...)
#4. 

