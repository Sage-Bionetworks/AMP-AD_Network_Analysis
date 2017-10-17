foo <- readLines('~/Desktop/amp_ad_top_targets')

source('enrichmentAnalysis/run_amp_ad_enrichment.R')

testList <- list(ad=foo,dummy=c('APOE','VGF'))

run_test <- run_amp_ad_enrichment(testList,'andhow',manifestId='syn11182793')
run_test <- dplyr::filter(run_test,category=='ad')

set.seed(2)
run_test$fisherPval[run_test$fisherPval==1] <- runif(length(which(run_test$fisherPval==1)),0,1)

run_full <- run_amp_ad_enrichment(testList,'whaaa')


run_full <- dplyr::filter(run_full,category=='ad')
run_full$fisherPval[run_full$fisherPval==1] <- runif(length(which(run_full$fisherPval==1)),0,1)


