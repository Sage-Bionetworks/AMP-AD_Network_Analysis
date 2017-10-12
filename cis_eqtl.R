synapseClient::synapseLogin()
cis_eqtl <- function(ensg='ENSG00000134516',synId = 'syn9726160'){
  ###query table for gene
  ###identify most significant snp
 # pb$tick()
  queryString <- paste0('select * from ',synId,' where ensemblId = \'', ensg,'\'')
  foobar <- synapseClient::synTableQuery(queryString)@values
  
  res <- c('chromosome' = NA,
           'location' = NA,
           'indexSNP'=NA,
           'distanceToTSS'=NA,
           'expressionIncreasingAllele'=NA,
           'qvalue' = NA,
           'gene' = NA,
           'ensembl' = ensg,
           'effect' = NA)
  
  if (nrow(foobar)>0){
    foobar <- dplyr::arrange(foobar,qvalue)
    if(foobar$qvalue[1]<0.05){
      res <- c('chromosome' = foobar$chromosome[1],
               'location' = foobar$location[1],
               'indexSNP'=foobar$snpId[1],
               'distanceToTSS'=foobar$distanceToTSS[1],
               'expressionIncreasingAllele'=foobar$expressionIncreasingAllele[1],
               'qvalue'=foobar$qvalue[1],
               'gene' = foobar$hugoId[1],
               'ensembl' = foobar$ensemblId[1],
               'effect' = foobar$effectEstimate[1])
      print(res)
    } 
  }
  return(res)

}
load('aggregate_module_mainfest.rda')
# getTopHobs <- function(x){
#   return(dplyr::arrange(x$anno,desc(hubs))$GeneID[1:10])
# }
# 
# topHubs <- sapply(fullManifest,getTopHobs)
# topHubs2 <- sapply(fullManifest,getTopHobs)
# topHubs<-names(sort(table(c(topHubs)),decreasing=T))
# 
# 
# 
# 
# blue <- fullManifest$aggregatePHGturquoisePHG$anno$GeneID

all_genes <- sapply(fullManifest,function(x){return(x$anno$GeneID)})
all_genes <- unlist(all_genes)
all_genes <- unique(all_genes)

no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)
parallel::clusterCall(cl,synapseClient::synapseLogin)
#fob <- names(sort(igraph::degree(foo),decreasing=T))[1:20]
#pb <- progress::progress_bar$new(total=length(blue))
foo1<- parallel::parSapply(cl,all_genes,cis_eqtl)
foo1 <- t(foo1)
foo1 <- data.frame(foo1,stringsAsFactors=F)
foo1 <- na.omit(foo1)
save(foo1,file='eqtl_res.rda')
igap_res <- rSynapseUtilities::loadDelimIntoList('syn10008574')
igap_res <- igap_res[[1]]
foo2 <- dplyr::filter(igap_res,
              MarkerName%in%foo1$indexSNP)

foo4 <- dplyr::left_join(foo2,foo1,by=c('MarkerName'='indexSNP'))




fxn1 <- function(x,foo4,igap_res){
  return(ks.test(foo4$Pvalue[foo4$ensembl%in%x$anno$GeneID],igap_res$Pvalue,alternative='greater')$p.value)
}
modPvalues <- sapply(fullManifest,fxn1,foo4,igap_res)
ks.test(foo4$Pvalue[foo4$ensembl%in%fullManifest$aggregateCBEblueCBE$anno$GeneID],igap_res$Pvalue,alternative='greater')

dim(foo4)


source('enrichmentAnalysis/run_amp_ad_enrichment.R')

gl <- list(cisqtl = foo1$ensembl,dummy = c('ENSG00000063322'))
foobar <- run_amp_ad_enrichment(gl,'ciseqtl',hgnc=F,manifestId='syn11182793')
View(foobar)


igap_res2 <- rSynapseUtilities::loadDelimIntoList('syn10008575')
igap_res2 <- igap_res2[[1]]
foo3 <- dplyr::filter(igap_res2,
                      MarkerName%in%foo1$indexSNP)


