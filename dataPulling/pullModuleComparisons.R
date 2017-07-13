synapseClient::synapseLogin()
foo <- synapseClient::synTableQuery("SELECT * FROM syn10163906")
bar <- foo@values
View(bar)
dlpfc1 <- grep('DLPFC',bar$ModuleNameFull)
dlpfc2 <- grep('DLPFC',bar$category)
keep <- intersect(dlpfc1,dlpfc2)
bar <- bar[keep,]
moduleNetwork <- igraph::graph_from_data_frame(bar[,c('ModuleNameFull','category')],directed = F)
moduleNetwork <- igraph::simplify(moduleNetwork)
moduleNetworkDf <- igraph::as_data_frame(moduleNetwork)
colnames(moduleNetworkDf) <- c('ModuleNameFull','category')
moduleNetworkDf2 <- dplyr::left_join(moduleNetworkDf,dplyr::select(bar,ModuleNameFull,category,fisherOR))
write.csv(moduleNetworkDf2,file='~/Desktop/rosmap.csv',quote=F)


#create module manifest with ad enrichment and cell type enrichments

genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))

source('enrichmentAnalysis/run_amp_ad_enrichment.R')
system.time(adTest <- run_amp_ad_enrichment(adList2,
                                            'ad_gwas',
                                            manifestId='syn10163855'))
adTest <- dplyr::filter(adTest,category=='ad_gwas')
adTestWide <- dplyr::select(adTest,ModuleNameFull,category,fisherPval) %>% tidyr::spread(category,fisherPval)

cellTypeTest <- run_amp_ad_enrichment(GeneSets$Cell_Markers,'cell_markers',manifestId='syn10163855')
wideCellType <- dplyr::select(cellTypeTest,ModuleNameFull,category,fisherPval) %>% tidyr::spread(category,fisherPval)

degs <- synapseClient::synGet('syn10163527') %>% synapseClient::getFileLocation()
load(degs)
degsEnrich <- run_amp_ad_enrichment(geneSet.test,'degs',manifestId='syn10163855')
wideDegsEnrich <- dplyr::select(degsEnrich,ModuleNameFull,category,fisherPval) %>% tidyr::spread(category,fisherPval)

View(wideDegsEnrich)

###build manifest
moduleManifest <- synapseClient::synTableQuery("SELECT * FROM syn10163855")@values
moduleManifest <- dplyr::select(moduleManifest,ModuleNameFull,method,brainRegion)
moduleManifestU <- moduleManifest[which(!duplicated(moduleManifest)),]

moduleManifestU <- dplyr::left_join(moduleManifestU,adTestWide)
moduleManifestU <- dplyr::left_join(moduleManifestU,wideCellType)
moduleManifestU <- dplyr::left_join(moduleManifestU,wideDegsEnrich)
colnames(moduleManifestU) <- gsub('\\.','_',colnames(moduleManifestU))
rSynapseUtilities::makeTable(moduleManifestU,"enrichment analyses July 8 2017","syn2370594")
dlpfcmanifest <- dplyr::filter(moduleManifestU,brainRegion=='DLPFC')
createCellTypeFactor <- function(x){
  minpval <- apply(x,1,min)
  wminpval <- apply(x,1,which.min)
  minpval<-p.adjust(minpval,method='fdr')
  wna <- which(minpval>0.05)
  dummyvar<-colnames(x)[wminpval]
  dummyvar[wna] <- NA
  return(dummyvar)
}
cellfactor <- createCellTypeFactor(dlpfcmanifest[,5:9])
dlpfcmanifest <- dplyr::mutate(dlpfcmanifest,cellfactor =createCellTypeFactor(dlpfcmanifest[,5:9]))
write.csv(dlpfcmanifest,file='~/Desktop/moduleanno.csv',quote=F)
write.csv(dlpfcmanifest[,c('ModuleNameFull','cellfactor')],file='~/Desktop/dlpfccells.csv',quote=F)

moduleManifestU<-dplyr::mutate(moduleManifestU,gwas_enrichment = -log10(p.adjust(moduleManifestU$ad_gwas,method='fdr')))


png(file='~/Desktop/adgwasenrich.png',
    height=800,
    width=1200,
    res=120,
    pointsize = 30)
g <- ggplot2::ggplot(moduleManifestU, 
                     ggplot2::aes(x=brainRegion,
                                  y=gwas_enrichment,
                                  fill=method))
g <- g + ggplot2::geom_col(position='dodge')
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g
dev.off()

fooSummarize <- dplyr::group_by(moduleManifestU,brainRegion,method)
fooSummarize2 <- dplyr::group_by(moduleManifestU,brainRegion,method) %>%
  dplyr::summarise(numberOfModules = length(method))


