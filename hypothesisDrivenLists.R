#construct lists for Thanneer


synapseClient::synapseLogin()
source('moduleAnalysis/summaryManifestFunctions.R')

#make ad lists and put into a table
adGeneticsSummaryAgg2 <- getAdGenetics2(synId='syn10915669')

#make cell lists and put into a table
genesets1 <- synapseClient::synGet('syn5923958')
load(genesets1@filePath)

#cell marker sets
GeneSets$Cell_Markers

cellDf<-utilityFunctions::list2df(GeneSets$Cell_Markers)
adDf <- utilityFunctions::list2df(adGeneticsSummaryAgg2[[2]])
combinedDf <- rbind(cellDf,adDf)
genesToMap <- unique(combinedDf$value)
mapTable <- utilityFunctions::convertHgncToEnsembl(genesToMap)

combinedDf<-dplyr::left_join(combinedDf,mapTable,by=c('value'='external_gene_name'))
View(combinedDf)
colnames(combinedDf) <- c('GeneSet','external_gene','ensembl_gene_id')
View(combinedDf)

rSynapseUtilities::makeTable(combinedDf,'Hypothesis Specific Gene Sets','syn2370594')
permLink <- githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                                 ref = 'branch',
                                 refName = 'module_ranking',
                                 repositoryPath = 'hypothesisDrivenLists.R')