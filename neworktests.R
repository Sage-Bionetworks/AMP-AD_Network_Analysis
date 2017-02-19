require(synapseClient)
synapseLogin()

foo <- synGet('syn8276546')
load(foo@filePath)
plot(bicNetworks$bicPath)
net <- as.matrix(bicNetworks$network)
hubs <- rowSums(net)
foo <- data.frame(ensembl_gene_id=names(hubs),
                     nEdges=hubs,
                     stringsAsFactors = F)

ensemblToGeneId <- function(ensemblId,dataSet='hsapiens_gene_ensembl'){
  require(dplyr)
  require(biomaRt)
  ensembl=useMart("ensembl")
  ensembl = useDataset(dataSet,mart=ensembl)
  gene1 <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),filters='ensembl_gene_id',values=ensemblId,mart=ensembl)
  return(gene1)
}
bar <- ensemblToGeneId(foo$ensemble_gene_id)
library(dplyr)

foobar <- dplyr::left_join(foo,bar,"ensembl_gene_id")
foobar <- dplyr::arrange(foobar,desc(nEdges))
cat(foobar$external_gene_name[1:500],file='~/Desktop/wut.csv',sep='\n')


###rosmap
foo2 <- synQuery('select name,id from file where method==\'bic\' and projectId==\'syn2370594\'')
foo3 <- synGet('syn8268669')
load(foo3@filePath)
rosmapnet <- as.matrix(bicNetworks$network)
hubs <- rowSums(rosmapnet)
foo <- data.frame(ensembl_gene_id=names(hubs),
                  nEdges=hubs,
                  stringsAsFactors = F)
foobar2 <- dplyr::left_join(foo,foobar,by="ensembl_gene_id")
colnames(foobar2)[c(2,3)] <- c('nEdgesRosmap','nEdgesMayoTCX')
foobar2 <- dplyr::mutate(foobar2,avgEdges=(nEdgesRosmap)/2+(nEdgesMayoTCX)/2)
foobar2 <- dplyr::arrange(foobar2,desc(avgEdges))
cat(foobar2$external_gene_name[1:2000],file='~/Desktop/wut.csv',sep='\n')
