cat('logging into Synapse...\n')
synapseClient::synapseLogin()
#grab module definitions
cat('pulling modules...\n')
allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ","syn9770791"))@values

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

source('run_amp_ad_enrichment_ensg.R')

pairwiseMods <- run_amp_ad_enrichment(modulesLargeList,
                                      "pairwiseComparison")
