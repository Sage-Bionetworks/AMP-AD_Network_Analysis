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
library(synapser)
synapser::synLogin()

genesets1 <- synapser::synGet('syn5923958')


load(genesets1$path)
adTypeDf <- utilityFunctions::list2df(GeneSets$Alzheimers)
write.csv(adTypeDf,file='~/Desktop/adTypeAnnotations.csv',quote=F)
write.csv(dplyr::filter(adTypeDf,key=='MouseMicroglia:2month_TG_vs_WT'),file='~/Desktop/two_months.csv',quote=F,row.names=F)
write.csv(dplyr::filter(adTypeDf,key=='MouseMicroglia:4month_TG_vs_WT'),file='~/Desktop/four_months.csv',quote=F,row.names=F)
write.csv(dplyr::filter(adTypeDf,key=='MouseMicroglia:6month_TG_vs_WT'),file='~/Desktop/six_months.csv',quote=F,row.names=F)
write.csv(dplyr::filter(adTypeDf,key=='MouseMicroglia:8month_TG_vs_WT'),file='~/Desktop/eight_months.csv',quote=F,row.names=F)
write.csv(dplyr::filter(adTypeDf,key=='MouseMicroglia:4month_vs_2month_TG-WT'),file='~/Desktop/four_two_months.csv',quote=F,row.names=F)


########
#get rosmap microglial/endothelial module
rosmap_mods <- synapser::synTableQuery("SELECT * FROM syn10309369 WHERE ModuleNameFull = \'metanetwork5DLPFC\' or ModuleNameFull = \'metanetwork19DLPFC\' or ModuleNameFull = \'metanetwork9DLPFC\'")
View(rosmap_mods$asDataFrame())


metanetwork5 <- dplyr::filter(rosmap_mods$asDataFrame(),ModuleNameFull=='metanetwork5DLPFC')
metanetwork19 <- dplyr::filter(rosmap_mods$asDataFrame(),ModuleNameFull=='metanetwork19DLPFC')
metanetwork9 <- dplyr::filter(rosmap_mods$asDataFrame(),ModuleNameFull=='metanetwork9DLPFC')
View(metanetwork5)
View(metanetwork9)


mods <- synapser::synTableQuery("SELECT * FROM syn10309369 where method =\'metanetwork\' and brainRegion=\'DLPFC\'")

mods <- mods$asDataFrame()
mods <- dplyr::select(mods,GeneID,external_gene_name)

View(foobar)
foobar_5 <- dplyr::filter(foobar,from %in% metanetwork5$GeneID & to %in% metanetwork5$GeneID)
foobar_5
foobar_5 <- dplyr::left_join(foobar_5,mods,c('from'='GeneID'))
foobar_5 <- dplyr::left_join(foobar_5,mods,c('to'='GeneID'))
write.csv(foobar_5,file='~/Desktop/mn5.csv',quote=F,row.names=F)
View(foobar_5)

foobar_9 <- dplyr::filter(foobar,from %in% metanetwork9$GeneID & to %in% metanetwork9$GeneID)
foobar_9
foobar_9 <- dplyr::left_join(foobar_9,mods,c('from'='GeneID'))
foobar_9 <- dplyr::left_join(foobar_9,mods,c('to'='GeneID'))
write.csv(foobar_9,file='~/Desktop/mn9.csv',quote=F,row.names=F)
View(foobar_9)


foobar_19 <- dplyr::filter(foobar,from %in% metanetwork19$GeneID & to %in% metanetwork19$GeneID)
foobar_19
foobar_19 <- dplyr::left_join(foobar_19,mods,c('from'='GeneID'))
foobar_19 <- dplyr::left_join(foobar_19,mods,c('to'='GeneID'))
write.csv(foobar_19,file='~/Desktop/mn19.csv',quote=F,row.names=F)
View(foobar_19)


