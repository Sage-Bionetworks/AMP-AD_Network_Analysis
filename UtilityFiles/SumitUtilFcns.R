Combine.DF.InFolder <- function(FolderName){
  
  FileNames <- list.files(FolderName)
  
  tempName <- paste(c(FolderName,'/',FileNames[1]),collapse = '')
  
  load(tempName)
  DF <- GeneFrame

  for (i in 2:length(FileNames)){
    tempName <- paste(c(FolderName,'/',FileNames[i]),collapse = '')
    load(tempName)
    DF <- rbind(DF, GeneFrame)
  }
  
  return(DF)
  
} 