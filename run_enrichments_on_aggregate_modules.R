synapseClient::synapseLogin()


#load enrichment functions

source('enrichmentAnalysis/run_amp_ad_enrichment.R')


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


source('enrichmentAnalysis/run_amp_ad_enrichment.R')
enrichr <- list()
enrichr$go_bp <- run_amp_ad_enrichment(go_bp,'go_bp',manifestId = 'syn11182793')
enrichr$go_cc <- run_amp_ad_enrichment(go_cc,
                                       'go_cc',
                                       manifestId = 'syn11182793')
enrichr$go_mf <- run_amp_ad_enrichment(go_mf,
                                       'go_mf',
                                       manifestId = 'syn11182793')
enrichr$kegg <- run_amp_ad_enrichment(kegg,
                                      'kegg',
                                      manifestId = 'syn11182793')
enrichr$wikipathways <- run_amp_ad_enrichment(wikipathways,
                                              'wikipathways',
                                              manifestId = 'syn11182793')
enrichr$reactome <- run_amp_ad_enrichment(reactome,
                                          'reactome',
                                          manifestId = 'syn11182793')
enrichr$biocarta <- run_amp_ad_enrichment(biocarta,
                                          'biocarta',
                                          manifestId = 'syn11182793')
enrichr$nci_nature <- run_amp_ad_enrichment(nci_nature,
                                            'nci_nature',
                                            manifestId = 'syn11182793')
enrichr$panther <- run_amp_ad_enrichment(panther,
                                         'panther',
                                         manifestId = 'syn11182793')
enrichr$ppi_hub_proteins <- run_amp_ad_enrichment(ppi_hub_proteins,'ppi_hub_proteins',manifestId = 'syn11182793')
enrichr$human_phenotype_ontology <- run_amp_ad_enrichment(human_phenotype_ontology,'human_phenotype_ontology',manifestId = 'syn11182793')
enrichr$jensen_diseases <- run_amp_ad_enrichment(jensen_diseases,'jensen_diseases',manifestId = 'syn11182793')
enrichr$omim_disease <- run_amp_ad_enrichment(omim_disease,'omim_disease',manifestId = 'syn11182793')
enrichr$omim_expanded <- run_amp_ad_enrichment(omim_expanded,'omim_expanded',manifestId = 'syn11182793')
enrichr$dbgap <- run_amp_ad_enrichment(dbgap,'dbgap',manifestId = 'syn11182793')
enrichr2 <- do.call(rbind,enrichr)

rSynapseUtilities::makeTable(enrichr2,'aggregate module pathway enrichments',projectId = 'syn2370594')
