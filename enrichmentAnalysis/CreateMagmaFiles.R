Create.MAGMA.GeneLists <- function(ModuleList, OutputFileName,
                                   manifestId = "syn9770791"){

  cat('pulling modules...\n')
  

  ModList <- list()
  for(i in 1:length(ModuleList))
  {
    cat('Loading module ',i,' of ',length(ModuleList),'\n')  
    temp <- synapseClient::synTableQuery(
      paste0("SELECT external_gene_name FROM ",manifestId," WHERE ModuleNameFull = ",
             "\'",ModuleList[i],"\'", sep= ""))@values
    
    #cat('List size',length(temp$external_gene_name),'\n')
    ModList <- append(ModList, as.list(ModuleList[i]))
    ModList <- as.list(append(ModList,
                              temp$external_gene_name))
    ModList <- append(ModList, as.list('\n'))
    #typeof(ModList)
    #Temp <- do.call(cbind, ModList)
    #cat(Temp, sep = " ")
    #cat("\n")
    
    
  }
  
  #sink()
  #writeLines(ModList,OutputFileName)
  #sink(OutputFileName)
  #print(ModList)
  #sink()
  #closeAllConnections()
  Temp <- paste(ModList, collapse = ' ')
  writeLines(Temp, OutputFileName)
  
  
  
  return(NULL)
}