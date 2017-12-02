#PLOT FOR HONG
#get bic dlpfc network
#get weights

#bic network: syn8268669
#rank consensus: syn8268680

synapseClient::synapseLogin()
foobar <- synapseClient::synGet('syn8268669')
load(foobar@filePath)
foobar3 <- synapseClient::synGet('syn8268680')
rankMatrix <- data.table::fread(foobar3@filePath,data.table=F)
rownames(rankMatrix) <- rankMatrix$V1
rankMatrix <- rankMatrix[,-1]
gc()
newMatrix <- as.matrix(rankMatrix)*bicNetworks$network
dim(newMatrix)
foo <- igraph::graph_from_adjacency_matrix(data.matrix(newMatrix),mode='undirected',weighted=T) %>%
  igraph::as_data_frame()
write.csv(foo,file='temp.csv',quote=F,row.names = F)

bar <- synapseClient::File('temp.csv',parentId='syn4213454')
bar <- synapseClient::synStore(bar)


synapser::synLogin()
bar <- synapser::synGet('syn11528498')

foobar <- data.table::fread(bar$path,data.table=F)
dim(foobar)
foobar$weight <- foobar$weight^(nrow(foobar)/4)
hist(foobar$weight)
sum(foobar$weight>.1)
foobar <- dplyr::filter(foobar,weight>.1)
mappingDf <- utilityFunctions::convertEnsemblToHgnc(unique(c(foobar$from,foobar$to)))
foobar<-dplyr::left_join(foobar,mappingDf,by=c('from'='ensembl_gene_id'))
foobar<-dplyr::left_join(foobar,mappingDf,by=c('to'='ensembl_gene_id'))

write.csv(foobar,file='~/Desktop/dlpfc.csv',quote=F)

#write lilly deg results to table
synapser::synLogin()

genesets1 <- synapser::synGet('syn5923958')


load(genesets1$path)
adTypeDf <- utilityFunctions::list2df(GeneSets$Alzheimers)
write.csv(adTypeDf,file='adTypeAnnotations.csv',quote=F)

