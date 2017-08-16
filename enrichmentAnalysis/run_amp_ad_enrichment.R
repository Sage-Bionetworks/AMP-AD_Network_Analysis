#########run enrichment function for manifest
run_amp_ad_enrichment <- function(geneSetList,
                                  geneSetName,
                                  hgnc = TRUE,
                                  manifestId = "syn10338156"){
  #INPUT:
  #geneSetList - a list of genes in hgnc or ensembl format
  #geneSetName - name of geneset (should be a single character string)
  #hgnc - boolean indicating whether the gene identifiers are hgnc (TRUE) or ensembl (FALSE)
  #manifestId - location of a synapse table with the modules definitions in terms of hgnc names


  #OUTPUT:
  #data frame with module name, gene set name, p-value, and odds ratio of enrichment
  library(dplyr)
  cat('logging into Synapse...\n')
  synapseClient::synapseLogin()
  #grab module definitions
  cat('pulling modules...\n')
  allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ",manifestId))@values

  listify <- function(x,y,z){
    ###fxn will listify a long form table
    ###x: unique key
    ###y: values
    ###z: keys
    return(unique(y[which(z==x)]))
  }
  cat('building module gene sets...\n')
  if(hgnc){
    modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                               listify,
                               allMods$external_gene_name,
                               allMods$ModuleNameFull)
  }else{
    modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                               listify,
                               allMods$GeneID,
                               allMods$ModuleNameFull)
  }
  names(modulesLargeList) <- unique(allMods$ModuleNameFull)
  cat('removing genes that are not relevant from reference set...\n')
  #get unique gene keys,drop categories in both cases that are 0 in size
  uniqueModuleList <- modulesLargeList %>%
    unlist %>%
    unique

  uniqueGeneSet <- geneSetList %>%
    unlist %>%
    unique

  refGeneSet <- uniqueModuleList

  cat('running enrichments....\n')

  res <- list()
  res$fisher <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                            modulesLargeList,
                                            geneSetList,
                                            refGeneSet)
  #pvalues are odd rows, odds ratios are even rows
  res$pval <- res$fisher[which(1:nrow(res$fisher)%%2==1),]
  rownames(res$pval) <- names(geneSetList)
  res$OR <- res$fisher[which(1:nrow(res$fisher)%%2==0),]
  rownames(res$OR) <- names(geneSetList)

  cat('producing tidy data frame....\n')
  #transpose pvalues
  pval1 <- t(res$pval)
  #make into a data frame
  pval1 <- data.frame(pval1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  pval1 <- dplyr::mutate(pval1,ModuleNameFull = rownames(pval1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest1 <- tidyr::gather(pval1,category,fisherPval,-ModuleNameFull)

  #transpose odds ratios
  or1 <- t(res$OR)
  #make into a data frame
  or1 <- data.frame(or1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  or1 <- dplyr::mutate(or1,ModuleNameFull = rownames(or1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest2 <- tidyr::gather(or1,category,fisherOR,-ModuleNameFull)

  #do a left join to combine the pvalues and odds ratios
  gatherTest <- dplyr::left_join(gatherTest1,
                                 gatherTest2)

  #add in the geneset names
  gatherTest <- dplyr::mutate(gatherTest,
                              geneSet = geneSetName)

  #return the data frame
  return(gatherTest)

}




run_amp_ad_enrichment_subset <- function(geneSetList,
                                  geneSetName,
                                  method,
                                  brainRegion,
                                  hgnc = TRUE,
                                  manifestId = "syn10338156"){

  #INPUT:
  #geneSetList - a list of genes in hgnc or ensembl format
  #geneSetName - name of geneset (should be a single character string)
  #hgnc - boolean indicating whether the gene identifiers are hgnc (TRUE) or ensembl (FALSE)
  #manifestId - location of a synapse table with the modules definitions in terms of hgnc names


  #OUTPUT:
  #data frame with module name, gene set name, p-value, and odds ratio of enrichment
  library(dplyr)
  cat('logging into Synapse...\n')
  #synapseClient::synapseLogin()
  #grab module definitions
  cat('pulling modules...\n')
  allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ",manifestId,
                                                " WHERE method = \'", method,
                                                "\' AND brainRegion = \'",
                                                brainRegion,"\'", sep=""))@values

  listify <- function(x,y,z){
    ###fxn will listify a long form table
    ###x: unique key
    ###y: values
    ###z: keys
    return(unique(y[which(z==x)]))
  }
  cat('building module gene sets...\n')
  if(hgnc){
    modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                               listify,
                               allMods$external_gene_name,
                               allMods$ModuleNameFull)
  }else{
    modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                               listify,
                               allMods$GeneID,
                               allMods$ModuleNameFull)
  }
  names(modulesLargeList) <- unique(allMods$ModuleNameFull)
  cat('removing genes that are not relevant from reference set...\n')
  #get unique gene keys,drop categories in both cases that are 0 in size
  uniqueModuleList <- modulesLargeList %>%
    unlist %>%
    unique

  uniqueGeneSet <- geneSetList %>%
    unlist %>%
    unique

  refGeneSet <- uniqueModuleList

  cat('running enrichments....\n')

  res <- list()
  res$fisher <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                                      modulesLargeList,
                                                      geneSetList,
                                                      refGeneSet)
  #pvalues are odd rows, odds ratios are even rows
  res$pval <- res$fisher[which(1:nrow(res$fisher)%%2==1),]
  rownames(res$pval) <- names(geneSetList)
  res$OR <- res$fisher[which(1:nrow(res$fisher)%%2==0),]
  rownames(res$OR) <- names(geneSetList)

  cat('producing tidy data frame....\n')
  #transpose pvalues
  pval1 <- t(res$pval)
  #make into a data frame
  pval1 <- data.frame(pval1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  pval1 <- dplyr::mutate(pval1,ModuleNameFull = rownames(pval1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest1 <- tidyr::gather(pval1,category,fisherPval,-ModuleNameFull)

  #transpose odds ratios
  or1 <- t(res$OR)
  #make into a data frame
  or1 <- data.frame(or1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  or1 <- dplyr::mutate(or1,ModuleNameFull = rownames(or1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest2 <- tidyr::gather(or1,category,fisherOR,-ModuleNameFull)

  #do a left join to combine the pvalues and odds ratios
  gatherTest <- dplyr::left_join(gatherTest1,
                                 gatherTest2)

  #add in the geneset names
  gatherTest <- dplyr::mutate(gatherTest,
                              geneSet = geneSetName)

  #return the data frame
  return(gatherTest)

}



Assign.Enriched.GS.Module <- function(DF){

  #get the unique module names in the dataframe
  ModNames <- unique(DF$ModuleNameFull)

  

  RetVar <- list()

  for (i in 1:length(ModNames)){

    #Find all indices for this module
    In <- which(DF$ModuleNameFull %in% ModNames[i])

    #Get the index of the minimum p-value
    In_min <- In[which.min(DF$fisherPval[In])]


    #Assign enrichment set to the module
    RetVar$EnrichCat[i] <- DF$category[In_min]

    #Assign minimum Pval
    RetVar$Pval[i] <- DF$fisherPval[In_min]

  }

  return(RetVar)

}
