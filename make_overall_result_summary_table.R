load('masterSheetEnrich.rda')
keep_cols <- c('AD_CONTROL_FEMALEDOWN',
               'AD_CONTROL_MALEDOWN',
               'AD_CONTROLDOWN',
               'AD_CONTROL_MALEUP',
               'AD_CONTROL_FEMALEUP',
               'AD_CONTROLUP',
               'Zhang.Endothelial',
               'Zhang.Astrocyte',
               'Zhang.MyelinOligos',
               'Zhang.Neuron',
               'Zhang.Microglia')

keep_rows <- c('aggregateDLPFCblackDLPFC',
               'aggregateDLPFCbrownDLPFC',
               'aggregatePHGgreenPHG',
               'aggregatePHGturquoisePHG',
               'aggregateDLPFCredDLPFC',
               'aggregateDLPFCyellowDLPFC',
               'aggregatePHGbluePHG',
               'aggregatePHGbrownPHG',
               'aggregatePHGyellowPHG',
               'aggregateCBEyellowCBE',
               'aggregateTCXblueTCX',
               'aggregateTCXturquoiseTCX',
               'aggregateCBEredCBE',
               'aggregateTCXredTCX',
               'aggregateCBEblueCBE',
               'aggregateCBEbrownCBE',
               'aggregateDLPFCturquoiseDLPFC',
               'aggregateTCXbrownTCX',
               'aggregateTCXpinkTCX',
               'aggregateTCXyellowTCX')
#keep_cols2 <- c('Zhang.Microglia',)
pheatmap::pheatmap((data.matrix(masterSheet[keep_rows,keep_cols])),cluster_rows = F,cluster_cols=F)
f3 <-(Matrix::Matrix(as.matrix(masterSheet[keep_rows,keep_cols])))

pheatmap::pheatmap(data.matrix(masterSheet[,keep_cols]))


module_cats <- list()
module_cats$female_up_no_male<-masterSheet$AD_CONTROL_FEMALEUP & !(masterSheet$AD_CONTROL_MALEDOWN | masterSheet$AD_CONTROL_MALEUP)
module_cats$female_down_no_male <- masterSheet$AD_CONTROL_FEMALEDOWN & !(masterSheet$AD_CONTROL_MALEDOWN | masterSheet$AD_CONTROL_MALEUP)
module_cats$female_down_male_down <- masterSheet$AD_CONTROL_FEMALEDOWN & masterSheet$AD_CONTROL_MALEDOWN
module_cats$female_up_male_up <- masterSheet$AD_CONTROL_FEMALEUP & masterSheet$AD_CONTROL_MALEUP
#module_cats$updownupdown <- masterSheet$AD_CONTROL_FEMALEUP & masterSheet$AD_CONTROL_FEMALEDOWN & masterSheet$AD_CONTROL_MALEUP & masterSheet$AD_CONTROL_MALEDOWN
ab<-unlist(sapply(module_cats,function(x,y){return(y[x])},rownames(masterSheet)))
f3 <-(Matrix::Matrix(as.matrix(masterSheet[ab,keep_cols])))
pheatmap::pheatmap(data.matrix(masterSheet[ab,keep_cols]),cluster_rows = F,cluster_cols = F)
#module_cats$updown <-  masterSheet$AD_CONTROL_FEMALEUP & masterSheet$AD_CONTROL_FEMALEDOWN
#module_cats$female_diff_male_diff <- masterSheet$AD_CONTROLDOWN & masterSheet$AD_CONTROLDOWN