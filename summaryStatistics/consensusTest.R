getTissueModules <- function(x = 'dorsolateralPrefrontalCortex'){
  string1 <- paste0("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'moduleIdentification' OR analysisType = 'consensusModuleIdentification') AND ( method = 'kmeans' OR method = 'wina' OR method = 'speakeasy' OR method = 'megena' OR method = 'wgcna' ) AND (tissueOfOrigin = '",x,"') AND (columnScaled = 'TRUE') )")
  
  moduleManifest <- synapseClient::synTableQuery(string1)
  
  return(moduleManifest@values)
}

synapseClient::synapseLogin()
foo <- getTissueModules()

### load delims into a list
bar <- rSynapseUtilities::loadDelimIntoList(foo$id)

### compute NMI and ARI