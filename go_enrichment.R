#gene ontology enrichment analyses

pullReferenceGeneSets <- function(url){
  #pull gene set from Enrichr
  foo <- RCurl::getURL(url)
  
  #split by carriage returns
  fooSplit <- strsplit(foo,"\n")
  
  #unlist
  fooSplit <- unlist(fooSplit)
  
  #split by tab
  fooSplit <- strsplit(fooSplit,"\t")
  
  #name them by the first element
  catNames <- lapply(fooSplit,function(x) x[1])
  names(fooSplit) <- catNames
  
  #drop the first two elements
  fooSplit <- lapply(fooSplit,function(x) x[-c(1,2)])
  
  return(fooSplit)
}

go_bp_2015_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2015"
go_bp_2015 <- pullReferenceGeneSets(go_bp_2015_url)

go_mf_2015_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2015"
go_mf_2015 <- pullReferenceGeneSets(go_mf_2015_url)

go_cc_2015_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2015"
go_cc_2015 <- pullReferenceGeneSets(go_cc_2015_url)


source('run_amp_ad_enrichment.R')
go <- list()
go$bp <- run_amp_ad_enrichment(go_bp_2015,
                               'go_bp_2015')
go$cc <- run_amp_ad_enrichment(go_cc_2015,
                               'go_cc_2015')
go$mf <- run_amp_ad_enrichment(go_mf_2015,
                               'go_mf_2015')
