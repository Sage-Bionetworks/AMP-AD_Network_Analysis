load('aggregate_module_mainfest.rda')

library(Matrix)
foobar<-utilityFunctions::convertAdjacencyToEdgeList(fullManifest$aggregatePHGturquoisePHG$adjacencyMatrix)
foobar1 <- utilityFunctions::convertEnsemblToHgnc(unique(c(foobar)))

foobar <- data.frame(foobar,stringsAsFactors=F)
View(foobar)
colnames(foobar) <- c('node1','node2')
foobar <- dplyr::left_join(foobar,foobar1,by=c('node1'='ensembl_gene_id'))
foobar <- dplyr::left_join(foobar,foobar1,by=c('node2'='ensembl_gene_id'))

write.csv(foobar[,3:4],file='phg_aggregate_turquoise.csv',quote=F)


coloc_genes <- data.frame(geneName = c('CR1','HLA-DRB5','PSMB9','HLA-DPB1','CD2AP','RIN3'),
                          value = rep(TRUE,6),
                          stringsAsFactors=F)
write.csv(coloc_genes,file='coloc_genes.csv',quote=F)
