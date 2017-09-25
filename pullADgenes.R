pullADgenes <- function(){
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  dbgap <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=dbGaP")
  #kegg, wikipathways, biocarta, panther, jensen_diseases, omim_disease, omim_expanded
  kegg <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016")
  wikipathways <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2016")
  jensen_diseases <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Jensen_DISEASES")
  biocarta <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioCarta_2016")
  panther <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")
  omim_disease <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Disease")
  omim_expanded <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Expanded")
  
  adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
  adList <- c(adList,'HLA-DRB5','HLA-DRB1')
  adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
  adList <- list(igap = adList)
  adList$dbgap <- dbgap$`Alzheimer Disease`
  adList$kegg <- kegg$`Alzheimer's disease_Homo sapiens_hsa05010`
  adList$wikipathwaysMouse <- wikipathways$`Alzheimers Disease_Mus musculus_WP2075`
  adList$wikipathwaysHuman <- wikipathways$`Alzheimers Disease_Homo sapiens_WP2059`
  adList$jensenDisease <- jensen_diseases$`Alzheimer's_disease`
  adList$biocarta <- biocarta$`Deregulation of CDK5 in Alzheimers Disease_Homo sapiens_h_p35alzheimersPathway`
  adList$pantherAmyloid <- panther$`Alzheimer disease-amyloid secretase pathway_Homo sapiens_P00003`
  adList$pantherPresenilin <- panther$`Alzheimer disease-presenilin pathway_Homo sapiens_P00004`
  adList$omim <- omim_disease$`alzheimer_disease`
  adList$omimExpanded <- omim_expanded$`alzheimer_disease`
  
  genecardsObj <- synapseClient::synGet('syn10507702')
  genecards <- data.table::fread(synapseClient::getFileLocation(genecardsObj),data.table=F)
  
  adList$genecards <- genecards$`Gene Symbol`
  return(adList)
}

pullDegGenes <- function(){
  degResObj <- synapseClient::synGet("syn10496554")
  load(degResObj@filePath)
  return(amp.ad.de.geneSets)
}

pullCellTypeGenes <- function(){
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  return(GeneSets$Cell_Markers)
}
