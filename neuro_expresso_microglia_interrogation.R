#load neuroexpresso gene sets
download.file("https://github.com/oganm/neuroExpressoAnalysis/raw/master/data/mouseMarkerGenesCombined.rda",'neuroex.rda')
load('neuroex.rda')

all_cell_types <- mapply(function(list_of_genes,name_of_list){
  name_of_lists <- names(list_of_genes)
  name_of_lists <- paste0(name_of_list,'_',name_of_lists)
  names(list_of_genes) <- name_of_lists
  return(list_of_genes)
},mouseMarkerGenesCombined,
names(mouseMarkerGenesCombined),
SIMPLIFY=T)

all_cell_types<-Reduce('c',all_cell_types)


#cortex_mg_act <- mouseMarkerGenesCombined$Cortex$Microglia_activation
#cortex_mg_deact <- mouseMarkerGenesCombined$Cortex$Microglia_deactivation
#pull aggregate modules

aggMods <- rSynapseUtilities::loadFullTable("syn10915669")

mouseIds <- utilityFunctions::getMouseOrthologsFromHumanIds(unique(aggMods$GeneID))

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


library(dplyr)
test_res <- utilityFunctions::fisherWrapperPval %>%
  utilityFunctions::outerSapplyParallel(dfList,
                                        all_cell_types,
                                        unique(aggMods$mmusculus_homolog_associated_gene_name))
test_res2 <- test_res*(nrow(test_res)*ncol(test_res))
test_res3<-(apply(test_res2<0.05,1,as.numeric))
cs_res <- colSums(test_res3)
rs_res <- rowSums(test_res3)
rownames(test_res3) <- colnames(test_res2)
test_res3 <- test_res3[which(rs_res>0),which(cs_res>0)]
pheatmap::pheatmap(t(test_res3))


