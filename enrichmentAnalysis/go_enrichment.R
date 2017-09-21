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
  
  fxn1 <- function(x){
    library(dplyr)
    strsplit(x,',') %>%
      sapply(function(y) y[1]) %>%
      return
  }
  #remove 1.0 from elements
  fooSplit2 <- lapply(fooSplit,fxn1)
  names(fooSplit2) <- names(fooSplit)
  return(fooSplit2)
}

go_bp_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2017b"
go_bp <- pullReferenceGeneSets(go_bp_url)
go_mf_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2017b"
go_mf <- pullReferenceGeneSets(go_mf_url)
go_cc_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2017b"
go_cc <- pullReferenceGeneSets(go_cc_url)
kegg_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016"
kegg <- pullReferenceGeneSets(kegg_url)
wikipathways <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2016")
reactome <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016")
biocarta <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioCarta_2016")
nci_nature <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=NCI-Nature_2016")
panther <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")
ppi_hub_proteins <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=PPI_Hub_Proteins")
human_phenotype_ontology <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Human_Phenotype_Ontology")
jensen_diseases <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Jensen_DISEASES")
omim_disease <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Disease")
omim_expanded <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Expanded")
dbgap <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=dbGaP")


source('run_amp_ad_enrichment.R')
enrichr <- list()
enrichr$go_bp <- run_amp_ad_enrichment(go_bp,
                               'go_bp')
enrichr$go_cc <- run_amp_ad_enrichment(go_cc,
                               'go_cc')
enrichr$go_mf <- run_amp_ad_enrichment(go_mf,
                               'go_mf')
enrichr$kegg <- run_amp_ad_enrichment(kegg,
                                      'kegg')
enrichr$wikipathways <- run_amp_ad_enrichment(wikipathways,
                                              'wikipathways')
enrichr$reactome <- run_amp_ad_enrichment(reactome,
                                          'reactome')
enrichr$biocarta <- run_amp_ad_enrichment(biocarta,
                                          'biocarta')
enrichr$nci_nature <- run_amp_ad_enrichment(nci_nature,
                                            'nci_nature')
enrichr$panther <- run_amp_ad_enrichment(panther,
                                         'panther')
enrichr$ppi_hub_proteins <- run_amp_ad_enrichment(ppi_hub_proteins,'ppi_hub_proteins')
enrichr$human_phenotype_ontology <- run_amp_ad_enrichment(human_phenotype_ontology,'human_phenotype_ontology')
enrichr$jensen_diseases <- run_amp_ad_enrichment(jensen_diseases,'jensen_diseases')
enrichr$omim_disease <- run_amp_ad_enrichment(omim_disease,'omim_disease')
enrichr$omim_expanded <- run_amp_ad_enrichment(omim_expanded,'omim_expanded')
enrichr$dbgap <- run_amp_ad_enrichment(dbgap,'dbgap')
