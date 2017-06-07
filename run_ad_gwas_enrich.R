synapseClient::synapseLogin()
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))

source('run_amp_ad_enrichment.R')
system.time(adTest <- run_amp_ad_enrichment(adList2,
                               'ad_gwas'))



system.time(aaaa <- utilityFunctions::outerSapply(utilityFunctions::fisherWrapper,
                                      modulesLargeList[1:100],
                                      modulesLargeList[1:30],
                                      unique(unlist(modulesLargeList))))

system.time(aaab <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                                  modulesLargeList[1:100],
                                                  modulesLargeList[1:30],
                                                  unique(unlist(modulesLargeList))))
parseOuterList <- function(x,ind){
  
}