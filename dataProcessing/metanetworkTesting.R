require(metanetwork)
require(synapseClient)
synapseLogin()
#get mayo network

#get rosmap network

foo <- synQuery('select name,id,study,tissueOfOrigin from file where projectId==\'syn2370594\' and method==\'bic\'')


bar <- sapply(foo$file.id,synGet)
names(bar) <- foo$file.tissueOfOrigin
load(bar$dorsolateralPrefrontalCortex@filePath)
rosmap <- as.matrix(bicNetworks$network)
load(bar$temporalCortex@filePath)
mayoTCX <- as.matrix(bicNetworks$network)
load(bar$cerebellum@filePath)
mayoCER <- as.matrix(bicNetworks$network)
dim(mayoCER)
dim(mayoTCX)
dim(rosmap)


######compare mayo tcx and rosmap
network1 <- rosmap
network2 <- mayoTCX
network3 <- mayoCER
geneName1 <- colnames(network1)
geneName2 <- colnames(network2)
geneName3 <- colnames(network3)
geneName <- intersect(geneName1,geneName2)
geneName <- intersect(geneName,geneName3)
network1 <- network1[geneName,geneName]
network2 <- network2[geneName,geneName]
network3 <- network3[geneName,geneName]

net1 <- network1[which(upper.tri(network1))]
net2 <- network2[which(upper.tri(network2))]
net3 <- network3[which(upper.tri(network3))]


foobar <- network1&network2
edgeList <- which(foobar,T)
edgeList2 <- edgeList
edgeList2[,1] <- colnames(foobar)[edgeList[,1]]
edgeList2[,2] <- colnames(foobar)[edgeList[,2]]

fb3 <- utilityFunctions::convertEnsemblToHgnc(unique(c(edgeList2[,1],edgeList2[,2])))
egn <- fb3$external_gene_name
edgeList3 <- edgeList2
edgeList3[,1] <- egn[edgeList2[,1]]
edgeList3[,2] <- egn[edgeList2[,2]]
write.csv(edgeList3,file='~/Desktop/edgeList.csv',quote=F)
names(egn) <- fb3$ensembl_gene_id
rosmapVmayoTCX <- table(net1,net2)
bar3 <- fisher.test(rosmapVmayoTCX)

library(VennDiagram)
venn.plot <- draw.triple.venn(sum(net1),sum(net2), sum(net3), sum(net1&net2),sum(net2&net3),sum(net1&net3),sum(net1&net2&net3), c("ROSMAP","Mayo TCX", "Mayo CER"))
grid.text("Edge Set Overlaps. OR=464,p-value<1e-16",.5,.9)
#grid.draw(venn.plot)

library(data.table)
tableStats <- fread('~/Desktop/tableStats.csv',data.table=F)
tableStats <- dplyr::select(tableStats,-SUID,-NumberOfDirectedEdges,-selected,-`shared name`)
tableStats <- dplyr::select(tableStats,-SelfLoops,-NumberOfUndirectedEdges,-PartnerOfMultiEdgedNodePairs,-IsSingleNode)

tableStats2 <- dplyr::select(tableStats,-name)%>%data.matrix
svd2 <- svd(scale(tableStats2))
tableStats <- dplyr::mutate(tableStats,scaledDegree = scale(as.numeric(Degree)))

tableStats$scaledStress <- scale(as.numeric(tableStats$Stress))
tableStats$degreeStress=tableStats$scaledStress/2+tableStats$scaledDegree/2
#tableStats <- dplyr::mutate(tableStats,scaledStress = scale(as.numeric(Stress)))


genesets1 <- synGet('syn5923958')
load(genesets1)
adList <- GeneSets$Alzheimers$`MouseMicroglia:2month_TG_vs_WT`
#adList <- scan('~/Desktop/adgenes.csv',what='character',sep=',')

#adList <- c(adList,'HLA-DRB5','HLA-DRB1')
alzGene <- tableStats$name%in%adList
library(glmnet)
foobar <- cv.glmnet(y=alzGene,x=scale(tableStats2),family='binomial')

summary(glm(alzGene~tableStats2[,c('Degree','NeighborhoodConnectivity')],family='binomial'))

