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

GWAS.Enrich.Modules <- function(ModName, AnnotFile, GWAS_file,
  manifestId = "syn9770791"){

  #Getting genes in the module
  cat('Getting genes in module...\n')
  #cat("SELECT external_gene_name FROM ",manifestId," WHERE ModuleNameFull = ",
  #       "\'",ModName,"\'", sep= "")
  temp <- synapseClient::synTableQuery(
    paste0("SELECT external_gene_name FROM ",manifestId," WHERE ModuleNameFull = ",
           "\'",ModName[1],"\'", sep= ""))@values

  GeneList <- temp$external_gene_name

  #reading the list of genes from the annotation file
  cat('Reading the list of genes from the annotation file ...\n')
  #Temp <- read.table(AnnotFile, skip = 2, header = FALSE, fill = TRUE)
  #Temp <- Temp$V1
  con <- file(AnnotFile, open = 'r')

  Temp <- c()

  while(length(oneLine <- readLines(con, n = 1, warn = F))>0){
    tmp <- unlist(strsplit(oneLine, split = "\t"))
    Temp <- c(Temp, tmp[1])
  }

  close(con)

  Temp <- Temp[3:length(Temp)]

  #reading the GWAS data
  cat('Reading the GWAS data ..\n')
  GWAS_data <- read.table(GWAS_file, header = TRUE)
  GWAS_snp <- GWAS_data$MarkerName
  GWAS_pval <- GWAS_data$Value

  rm(GWAS_data)

  Pval_list <- rep(1,length(GeneList))

  cat('Finding smallest p-value for each gene... \n')

  for ( i in 1:length(GeneList)){

    cat('Gene',i,'of',length(GeneList),'...\n')

    In <- which(Temp %in% GeneList[i])
    #print(In)

    if(length(In)==0){
      next
    }

    SNP_list <- read.table(AnnotFile, skip = In[1] + 2,
    header = FALSE, nrows = 1)

    cat('Number of SNPs is',length(SNP_list),'\n')

    SNP_list <- DF2List(SNP_list)

    #cat('Converted to list \n')

    #print(SNP_list)
    #Pval_In <- unlist(lapply(SNP_list, function(x) grep(x, GWAS_snp)))
    #print(typeof(Pval_In))
    #Pval_temp <- Pval_list[Pval_In]
    #print(Pval_temp)

    Pval_In <- Get.Min.PvalIn(SNP_list, GWAS_snp)
    #cat('Obtained index of minimum pval \n')
    Pval_temp <- GWAS_pval[Pval_In]

    Pval_list[i] <- min(Pval_temp)
    cat('Pval for gene',GeneList[i],'is',Pval_list[i],'\n')

  }

  RetVar <- list()
  RetVar$Genes <- GeneList
  RetVar$Pval_min <- Pval_list

  return(RetVar)

}

DF2List <- function(DF){
  #function to convert one line dataframe to vector
  DF.list <- as.list(DF)
  my.vec <- lapply(DF.list, function(x){x[[1]][[1]]})
  my.final.vec <- as.data.frame(unlist(my.vec), stringsAsfactors = F)
  FinalList <- as.vector(my.final.vec[,1])
  return(FinalList)
}

Get.Min.PvalIn <- function(SNP_list, GWAS_snp){
  #Pval_In <- unlist(lapply(SNP_list, function(x) grep(x, GWAS_snp)))
  #print(typeof(Pval_In))
  #print(Pval_In)
  Pval_In <- which(GWAS_snp %in% SNP_list)

  return(Pval_In)


}



Gen.GSEA.DF <- function(GSE_file, AnnotFile, GWAS_file){

  cat('Obtaining most enriched module in region .. \n')
  GWAS_enrich <- read.table(GSE_file, skip = 3, header = T)

  Min_In <- which.min(GWAS_enrich$P)
  Min_mod <- as.vector(droplevels(GWAS_enrich$SET[Min_In]))

  cat('Obtaining enriched genes in module .. \n')
  GenePvalList <- GWAS.Enrich.Modules(Min_mod,
    AnnotFile = AnnotFile, GWAS_file = GWAS_file)

  GenePvalList$MagmaPval <- rep(GWAS_enrich$P[Min_In],
    length(GenePvalList$Pval_min))
  GenePvalList$ModuleName <- rep(Min_mod,
    length(GenePvalList$Pval_min))
  cat('Returning output data frame ..\n')

  GeneFrame <- data.frame(GenePvalList)

  return(GeneFrame)

}
