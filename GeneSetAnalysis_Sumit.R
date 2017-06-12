#this function assumes that synapse has been logged into 
GetGenesFromModule <- function(ModuleName,
                               manifestId = "syn9770791"){
  
  library(dplyr)
  cat('logging into Synapse...\n')
  #grab module definitions
  cat('pulling genes...\n')
  cat("SELECT * FROM ",manifestId," WHERE ModuleNameFull = ",
      "\'",ModuleName,"\'", sep="")
  allMods <- synapseClient::synTableQuery(
    paste0("SELECT * FROM ",manifestId," WHERE ModuleNameFull = ",
           "\'",ModuleName,"\'", sep= ""))@values
  return(allMods)
}

GetModulesInClust <- function(ModuleList, ClustAssign, In){
  
  SubSetModule <- ModuleList[ClustAssign == In]
  return(SubSetModule)
  
}

GetUnionGenes <- function(ModuleList, manifestId = "syn9770791"){
  
  cat('pulling modules...\n')
  ModList <- list()
  for(i in length(ModuleList))
  {
    temp <- synapseClient::synTableQuery(
      paste0("SELECT GeneID FROM ",manifestId," WHERE ModuleNameFull = ",
             "\'",ModuleList[i],"\'", sep= ""))@values
  ModList <- append(ModList, temp)  
    
  }
  
  ModList <- unique(ModList)
  return(ModList)
}

