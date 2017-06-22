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


makeTable <- function(df,tableName,projectId){
  library(synapseClient)
  #synapseLogin()
  tcresult<-as.tableColumns(df)
  cols<-tcresult$tableColumns
  fileHandleId<-tcresult$fileHandleId
  project<- synGet(projectId)
  schema<-TableSchema(name=tableName, parent=project, columns=cols)
  table<-Table(schema, fileHandleId)
  table<-synStore(table, retrieveData=TRUE)
}



Create.GeneSet.List <- function(Gene_file){

  con <- file(Gene_file, open = 'r')

  GeneList <- list()

  cnt <- 1
  while(length(oneLine <- readLines(con, n = 1, warn = F))>0){
    tmp <- unlist(strsplit(oneLine, split = "\t"))
    tempStr <- paste(c('GeneList$mod',cnt,'<- tmp[3:length(tmp)]'),
    collapse = '')
    cnt <- cnt + 1
    #print(tempStr)
    #GeneList$tmp[1] <- tmp[3:length(tmp)]
    eval(parse(text = tempStr))
    #head(GeneList)
    #GeneList$cnt <- tmp[3:length(tmp)]

  }

  close(con)

  return(GeneList)


}
