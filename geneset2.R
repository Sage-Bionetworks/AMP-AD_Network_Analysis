#pull Bic network for ROSMAP
synapseClient::synapseLogin()
foo <- synapseClient::synGet('syn8268669')
load(foo@filePath)
bar <- which(bicNetworks$network!=0,T)
baz <- bar
baz[,1] <- colnames(bicNetworks$network)[bar[,1]]
baz[,2] <- colnames(bicNetworks$network)[bar[,2]]
write.csv(baz,file='~/Desktop/rosmap.csv',quote=F)
#pull gene set

metanetworkModules <- synapseClient::synTableQuery("SELECT * FROM syn10146524 WHERE ( ( method = 'metanetwork' ))")@values

synapseClient::synapseLogin()
geneSets1<-synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(geneSets1))
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
tl <- c(GeneSets$Alzheimers,GeneSets$Cell_Markers)

modulesLargeList <- lapply(unique(metanetworkModules$ModuleNameFull),
                           listify,
                           metanetworkModules$external_gene_name,
                           metanetworkModules$ModuleNameFull)
names(modulesLargeList) <- unique(metanetworkModules$ModuleNameFull)

aaaa <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,
                                                  modulesLargeList,
                                                  tl,
                                                  unique(unlist(modulesLargeList)))
bbbb <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,
                                      modulesLargeList,
                                      tl,
                                      unique(unlist(modulesLargeList)))
bbb <- bbbb < 0.05/(nrow(bbbb)*ncol(bbbb))
bbb <- apply(bbb,2,as.numeric)
rownames(bbb) <- rownames(bbbb)
aaaa <- aaaa[,which(colSums(bbb)>0)]
bbb <- bbb[,which(colSums(bbb)>0)]
aaaa <- aaaa*bbb
pheatmap::pheatmap(aaaa^.25)

#write bic network to cytoscape ready file

#import into cytoscape

#