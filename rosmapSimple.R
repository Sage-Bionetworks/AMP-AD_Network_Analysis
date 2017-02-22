library(synapseClient)
synapseLogin()

foo <- synGet('syn8268669')
load(foo@fileP)
library(utilityFunctions)
foo <- as.matrix(bicNetworks$network)
hubs <- rowSums(foo)
hubs2 <- data.frame(ensembl_gene_id=names(hubs),
                    nEdges=hubs,
                    stringsAsFactors=F)
bar <- utilityFunctions::convertEnsemblToHgnc(hubs2$ensembl_gene_id)
foobar <- dplyr::left_join(hubs2,bar,by='ensembl_gene_id')
foobar <- dplyr::arrange(foobar,desc(nEdges))
cat(foobar$external_gene_name[1:500],file='~/Desktop/testSet.csv',sep='\n')
foo3 <- which(foo!=0,T)
foo3 <- cbind(foo3,colnames(foo)[foo3[,'row']],colnames(foo)[foo3[,'col']])
colnames(foo3) <- c('n1','n2','ensembl_gene_id1','ensembl_gene_id2')
foo3 <- data.frame(foo3,stringsAsFactors=F)
foo3[1:5,]
foo3 <- merge(foo3,bar,by.x='ensembl_gene_id1',by.y='ensembl_gene_id')
colnames(foo3)[5] <- c('external_gene_name1')
foo3 <- merge(foo3,bar,by.x='ensembl_gene_id2',by.y='ensembl_gene_id')
edgeList <- dplyr::select(foo3,external_gene_name1,external_gene_name)
write.csv(edgeList,file='~/Desktop/rosmap.csv',quote=F,row.names=F)
maxEdges <- round((600*17000))/20
