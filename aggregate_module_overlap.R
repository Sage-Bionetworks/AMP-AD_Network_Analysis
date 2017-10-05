##login to synapse
synapseClient::synapseLogin()

##pull aggregate modules
aggregateModules <- rSynapseUtilities::loadFullTable('syn10915669')

View(aggregateModules)

mats <- utilityFunctions::pairwiseMatrixOfEnrichments(key=aggregateModules$ModuleNameFull,
                                                      value=aggregateModules$GeneID)

#get rid of underflow issues
mats$pval <- mats$pval + 1e-300

pheatmap::pheatmap(-log10(mats$pval))
#pheatmap::pheatmap(scale(t(scale(mats$or,center = F)),center=F))
