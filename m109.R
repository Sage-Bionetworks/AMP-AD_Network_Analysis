source('enrichmentAnalysis/run_amp_ad_enrichment.R')

synapseClient::synapseLogin()

foo <- synapseClient::synTableQuery("SELECT distinct geneName, hgncName FROM syn5321231 where speakeasyModule = \'109\'")@values


