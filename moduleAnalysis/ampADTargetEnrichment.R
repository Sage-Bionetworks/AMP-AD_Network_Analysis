synapseClient::synapseLogin()

foobar <- synapseClient::synGet('syn8656625')
foo <- data.table::fread(foobar@filePath,data.table=F)
source('enrichmentAnalysis/run_amp_ad_enrichment.R')
buh <- foo$ensembl_id
buh <- list()
buh$all <- foo$ensembl_id
buh$broad <- dplyr::filter(foo,group=='Broad-Rush')$ensembl_id
buh$duke <- dplyr::filter(foo,group=='Duke')$ensembl_id
buh$emory <- dplyr::filter(foo,group=='Emory')$ensembl_id
buh$mssm <- dplyr::filter(foo,group=='MSSM')$ensembl_id
buh$ufl <- dplyr::filter(foo,group=='UFL-ISB-Mayo')$ensembl_id
bar <- run_amp_ad_enrichment(buh,geneSetName = 'ampad',hgnc=FALSE,manifestId="syn10158502")
