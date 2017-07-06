cat('logging into Synapse...\n')
synapseClient::synapseLogin()
#grab module definitions
cat('pulling modules...\n')
allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ","syn10163855"))@values
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
cat('building module gene sets...\n')
modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                           listify,
                           allMods$GeneID,
                           allMods$ModuleNameFull)
names(modulesLargeList) <- unique(allMods$ModuleNameFull)

source('enrichmentAnalysis/run_amp_ad_enrichment.R')

pairwiseMods <- run_amp_ad_enrichment(modulesLargeList,
                                      "pairwiseComparison",
                                      hgnc=FALSE,
                                      manifestId = "syn10163855")
pairwiseMods <- dplyr::mutate(pairwiseMods,adj = p.adjust(fisherPval,method='fdr'))
pairwiseMods <- dplyr::filter(pairwiseMods,adj<=0.05)
rSynapseUtilities::makeTable(pairwiseMods,"pairwise module overlap manifest July 6 2017 FDR 0.05",'syn2370594')
