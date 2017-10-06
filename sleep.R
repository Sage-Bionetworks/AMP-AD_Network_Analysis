###being a morning person

##list derived from union of doi:10.1038/ncomms10889 and doi:10.1038/ncomms10448
morning_list <- list(sleep= c('RGS16','VIP','PER2','RASD1','PER3','FBXL3','PLCL1','APH1A','FBXL13','NOL4','TOX3','DLX5','RNASEL','VAMP3','CLN5','CA14','FAM185A','TNRC6B','RNF10','ADCY8','FAT1','ERC2','ASB1','AK5'),dummy=c('VGF'))


source('enrichmentAnalysis/run_amp_ad_enrichment.R')
foobar<-run_amp_ad_enrichment(morning_list,geneSetName = 'sleep',hgnc = T,manifestId = 'syn10338156')
foobar <- dplyr::filter(foobar,category=='sleep')
View(foobar)
