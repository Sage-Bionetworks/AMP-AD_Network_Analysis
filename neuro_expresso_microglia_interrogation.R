#load neuroexpresso gene sets
download.file("https://github.com/oganm/neuroExpressoAnalysis/raw/master/data/mouseMarkerGenesCombined.rda",'neuroex.rda')
load('neuroex.rda')

cortex_mg_act <- mouseMarkerGenesCombined$Cortex$Microglia_activation
cortex_mg_deact <- mouseMarkerGenesCombined$Cortex$Microglia_deactivation



library(utilityFunctions)
detach("package:utilityFunctions",force=T,unload=T)

#pull aggregate modules

aggMods <- rSynapseUtilities::loadFullTable("syn10915669")

mouseIds <- utilityFunctions::getMouseOrthologsFromHumanIds(unique(aggMods$GeneID))

View(mouseIds)


###run enrichments on aggregate modules for activated and deactivated gene sets

#make mouse version of module lists
aggMods <- dplyr::left_join(aggMods,
                            mouseIds,
                            by=c('GeneID' = "ensembl_gene_id"))

dfList <- lapply(unique(aggMods$ModuleNameFull),
                 utilityFunctions::listify,
                 aggMods$mmusculus_homolog_associated_gene_name,
                 aggMods$ModuleNameFull)
names(dfList) <- unique(aggMods$ModuleNameFull)

mg_act <- sapply(dfList,
                utilityFunctions::fisherWrapperPval,
                cortex_mg_act,
                unique(aggMods$mmusculus_homolog_associated_gene_name))

mg_deact <- sapply(dfList,
                             utilityFunctions::fisherWrapperPval,
                             cortex_mg_deact,
                             unique(aggMods$mmusculus_homolog_associated_gene_name))

spike_in <- mouseMarkerGenesCombined$Cortex$Astrocyte

test_res <- mg_act <- sapply(dfList,
                             utilityFunctions::fisherWrapperPval,
                             spike_in,
                             unique(aggMods$mmusculus_homolog_associated_gene_name))
