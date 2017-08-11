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

bbb2 <- bbb[-c(9,10),]

###get annotation columns
mani <- dplyr::select(metanetworkModules,Module,ModuleNameFull,brainRegion)
mani <- mani[!duplicated(mani),]
a2 <- mani$ModuleNameFull%in%colnames(bbb2)
mani2 <- mani[a2,]
brainRegion1 <- mani2$brainRegion
names(brainRegion1) <- mani2$ModuleNameFull
modulen <- mani2$Module
names(modulen) <- mani2$ModuleNameFull

brainRegion1 <- brainRegion1[colnames(bbb2)]
modulen <- modulen[colnames(bbb2)]
bbb3 <- bbb2
colnames(bbb3) <- modulen

df1 <- data.frame(brainRegion=brainRegion1,stringsAsFactors=F)
rownames(df1) <- names(brainRegion1)

png(file='~/Desktop/mousemganalysis.png',
    height=800,
    width=1200,
    res=120,
    pointsize = 20)
pheatmap::pheatmap(bbb2,color = c('white','blue'),border_color = NA,annotation_col=df1)
dev.off()
#write bic network to cytoscape ready file

#import into cytoscape

#