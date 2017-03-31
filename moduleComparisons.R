synapseClient::synapseLogin()

###get all label propagation modules

foo <- synapseClient::synQuery('select name,id,study from file where projectId==\'syn2370594\' and method==\'spinglass\' and analysisType==\'moduleIdentification\'')


####enrichment analyses for cell types
baz <- rSynapseUtilities::loadDelimIntoList(foo$file.id)
#sort(table(baz[[2]]$moduleNumber),decreasing=T)[1:20]

#add in hgnc data to modules
addHgnc <- function(df){
  geneTable <- utilityFunctions::convertEnsemblToHgnc(df$Gene.ID)
  df <- dplyr::left_join(df,geneTable,c('Gene.ID'='ensembl_gene_id'))
  return(df)
}

baz2 <- lapply(baz,addHgnc)
#convert to lists
convertToList <- function(df){
  tab1 <- table(df$moduleLabel)
  tab1 <- tab1[which(tab1>=20)]
  gn <- names(tab1)
  fxn1 <- function(x,df2){
    return(dplyr::filter(df2,moduleLabel==x)$external_gene_name)
  }
  listified <- lapply(gn,fxn1,df)
  names(listified) <- gn
  return(listified)
}

baz3 <- lapply(baz2,convertToList)

#get gene sets
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList


microglialList <- GeneSets$Cell_Markers$`Zhang:Microglia`
neuralList <- GeneSets$Cell_Markers$`Zhang:Neuron`
astrocyteList <- GeneSets$Cell_Markers$`Zhang:Astrocyte`
endothelialList <- GeneSets$Cell_Markers$`Zhang:Endothelial`

library(dplyr)





getPvaluesAndOddsRatios <- function(list1,list2,geneSet1){
  model <- list()
  model$msbbAdPval <- lapply(list1[[1]],
                       utilityFunctions::fisherWrapperPval,
                       geneSet1,
                       unique(baz2[[1]]$external_gene_name)) %>%
                       unlist %>%
                       p.adjust(method='fdr')
  
  model$msbbAdOR <- lapply(list1[[1]],
                     utilityFunctions::fisherWrapperOR,
                     geneSet1,
                     unique(list2[[1]]$external_gene_name)) %>%
                     unlist
  
  #rosmap
  model$rosmapAdPval <- lapply(list1[[2]],
                         utilityFunctions::fisherWrapperPval,
                         geneSet1,
                         unique(list2[[2]]$external_gene_name)) %>%
                         unlist %>%
                         p.adjust(method='fdr')
  
  model$rosmapAdOR <- lapply(list1[[2]],
                       utilityFunctions::fisherWrapperOR,
                       geneSet1,
                       unique(list2[[2]]$external_gene_name)) %>% 
                       unlist
  
  #mayo tcx
  model$mayoTcxAdPval <- lapply(list1[[3]],
                          utilityFunctions::fisherWrapperPval,
                          geneSet1,
                          unique(list2[[3]]$external_gene_name)) %>% 
                          unlist %>%
                          p.adjust(method='fdr')
  
  model$mayoTcxAdOR <- lapply(list1[[3]],
                        utilityFunctions::fisherWrapperOR,
                        geneSet1,
                        unique(list2[[3]]$external_gene_name)) %>% unlist
  
  return(model)
}

getPvaluesAndOddsRatios(baz3,baz2,microglialList)
#msbb



#utility function that pulls bic network, extracts things
rosmapNetworkObj <- synGet('syn8268669')
load(rosmapNetworkObj@filePath)
net1 <- as.matrix(bicNetworks$network)
tanRosmapManifest <- dplyr::filter(baz2[[2]],moduleLabel=='tan')
net1 <- net1[tanRosmapManifest$Gene.ID,tanRosmapManifest$Gene.ID]
net2 <- which(net1!=0,T)
net3 <- net2
net3[,1] <- tanRosmapManifest$external_gene_name[net2[,1]]
net3[,2] <- tanRosmapManifest$external_gene_name[net2[,2]]
write.csv(net3,file='~/Desktop/rosampAdGenetics.csv',quote=F)
cat(tanRosmapManifest$external_gene_name,file='~/Desktop/tan.csv',sep='\n')
foobar <- data.frame(external_gene_name=tanRosmapManifest$external_gene_name,
                     adLoci=tanRosmapManifest$external_gene_name%in%adList,
                     stringsAsFactors=F)
write.csv(foobar,file='~/Desktop/gene_names.csv',quote=F)
#run things
#utilityFunctions::fisherWrapperPval()
####comaprison of modules
rosmapVmsbbPval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,baz3[[2]],baz3[[1]],allGenes=union(unlist(baz3[[2]]),unlist(baz3[[1]])))
rosmapVmsbbOR <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,baz3[[2]],baz3[[1]],allGenes=union(unlist(baz3[[2]]),unlist(baz3[[1]])))

png(filename='~/Desktop/rosmap_msbb.png',height=4,width=6,units='in',pointsize = 5,res=300)
pheatmap::pheatmap(t((rosmapVmsbbOR)^(1/2)),
                   scale='none',
                   labels_col=names(baz3[[1]]),
                   xlab='ROSMAP DLPFC Modules',
                   ylab='MSBB FP Modules')
dev.off()

rosmapVmayoPval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,baz3[[2]],baz3[[3]],allGenes=union(unlist(baz3[[2]]),unlist(baz3[[3]])))
rosmapVmayoOR <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,baz3[[2]],baz3[[3]],allGenes=union(unlist(baz3[[2]]),unlist(baz3[[3]])))

png(filename='~/Desktop/rosmap_mayo.png',height=4,width=6,units='in',pointsize = 5,res=300)
pheatmap::pheatmap(t((rosmapVmayoOR)^(1/2)),
                   scale='none',
                   labels_col=names(baz3[[1]]),
                   xlab='ROSMAP DLPFC Modules',
                   ylab='Mayo TCX Modules')
dev.off()


msbbVmayoPval <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperPval,baz3[[1]],baz3[[3]],allGenes=union(unlist(baz3[[1]]),unlist(baz3[[3]])))
msbbVmayoOR <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapperOR,baz3[[1]],baz3[[3]],allGenes=union(unlist(baz3[[1]]),unlist(baz3[[3]])))

png(filename='~/Desktop/msbb_mayo.png',height=4,width=6,units='in',pointsize = 5,res=300)
pheatmap::pheatmap(t((msbbVmayoOR)^(1/2)),
                   scale='none',
                   labels_col=names(baz3[[1]]),
                   xlab='MSBB FP Modules',
                   ylab='Mayo TCX Modules')
dev.off()

#pheatmap::pheatmap(-log10(rosmapVmsbbPval+1e-200),scale = 'none')
####phenotype associations
w1 <- which(baz2[[2]]$external_gene_name%in%baz3[[2]]$tan)
ensgIdsRosmap <- baz2[[2]]$Gene.ID[w1]

combinedData[1:5,1:5]



microglia <- list(`ROSMAP\nDLPFC Tan`=baz3[[2]]$tan,`Mayo TCX\nLight Yellow`=baz3[[3]]$lightyellow,`MSSM FP Red`=baz3[[1]]$red)
VennDiagram::venn.diagram(microglia,filename='~/Desktop/test.png',sub.cex=.5,cat.cex=1,cat.dist=.18,height=3000,width=3000,resolution=800,imagetype="png",margin=0.05)


####enrichment for gwas associations