source('enrichmentAnalysis/run_amp_ad_enrichment.R')

synapseClient::synapseLogin()

foo <- synapseClient::synTableQuery("SELECT distinct geneName, hgncName FROM syn5321231 where speakeasyModule = \'109\'")@values

testSet <- list(m109=foo$hgncName, dummy=c('VGF','APOE'))
a1<-run_amp_ad_enrichment(testSet,'m109',manifestId='syn11182793')
View(a1)
