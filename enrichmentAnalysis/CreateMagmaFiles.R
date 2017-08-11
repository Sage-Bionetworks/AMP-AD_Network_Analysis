Create.MAGMA.GeneLists <- function(ModuleList, OutputFileName,
                                   manifestId = "syn9770791"){

  #The purpose of this function is to create a text file
  #containing the genes in a set of modules

  #Inputs:
  #ModuleList: A list/vector containing the module names
  #OutputFileName: The name of the output text file
  #manifestId: manifestId in synapse

  #Outputs:
  #None

  #pulling modules
  cat('pulling modules...\n')


  ModList <- list()
  for(i in 1:length(ModuleList))
  {
    #loading genes in each module from synapse
    cat('Loading module ',i,' of ',length(ModuleList),'\n')
    temp <- synapseClient::synTableQuery(
      paste0("SELECT external_gene_name FROM ",manifestId," WHERE ModuleNameFull = ",
             "\'",ModuleList[i],"\'", sep= ""))@values


    #appending genes in the extracted module to a list
    ModList <- append(ModList, as.list(ModuleList[i]))
    ModList <- as.list(append(ModList,
                              temp$external_gene_name))
    ModList <- append(ModList, as.list('\n'))


  }

  #writing the genes by module into a file
  Temp <- paste(ModList, collapse = ' ')
  writeLines(Temp, OutputFileName)

  return(NULL)
}







GWAS.Enrich.Modules <- function(ModName, AnnotFile, GWAS_file,
  manifestId = "syn9770791"){

  #The purpose of this function is to get the minimum SNP p-value
  #for all the genes in a module

  #Inputs:
  #ModName: Module name
  #AnnotFile: File containing the gene annotations for the GWAS_file
  #GWAS_file: File containing SNP p-values
  #manifestId: manifestId for synapse

  #Outputs:
  #RetVar$GeneList: The list of genes in the module
  #RetVar$Pval_min: The minimum p-value for each gene in the list


  #Getting genes in the module
  cat('Getting genes in module...\n')
  temp <- synapseClient::synTableQuery(
    paste0("SELECT external_gene_name FROM ",manifestId," WHERE ModuleNameFull = ",
           "\'",ModName[1],"\'", sep= ""))@values

  GeneList <- temp$external_gene_name

  #reading the list of genes from the annotation file
  cat('Reading the list of genes from the annotation file ...\n')
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

  #Finding the smallest p-value for each gene
  cat('Finding smallest p-value for each gene... \n')

  for ( i in 1:length(GeneList)){

    cat('Gene',i,'of',length(GeneList),'...\n')

    In <- which(Temp %in% GeneList[i])

    if(length(In)==0){
      next
    }

    SNP_list <- read.table(AnnotFile, skip = In[1] + 2,
    header = FALSE, nrows = 1)

    cat('Number of SNPs is',length(SNP_list),'\n')

    #Converting the dataframe object to list
    SNP_list <- DF2List(SNP_list)

    #Getting the minimum p-value from a given set of SNPs
    Pval_In <- Get.Min.PvalIn(SNP_list, GWAS_snp)
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

  #Inputs:
  #DF: A one line dataframe

  #Outputs:
  #FinalList: A list containing the same data

  DF.list <- as.list(DF)
  my.vec <- lapply(DF.list, function(x){x[[1]][[1]]})
  my.final.vec <- as.data.frame(unlist(my.vec), stringsAsfactors = F)
  FinalList <- as.vector(my.final.vec[,1])
  return(FinalList)
}







Get.Min.PvalIn <- function(SNP_list, GWAS_snp){

  #function to get the minimum p-value index for a given set of
  #SNPs from a larger document containing all the SNPs and their
  #p-values

  #Inputs:
  #SNP_list: A string list containing a subset of the total SNPs
  #GWAS_snp: The list of total SNPs along with their p-values

  #Outputs:
  #Pval_In: The index of the smallest SNP

  Pval_In <- which(GWAS_snp %in% SNP_list)

  return(Pval_In)

}



Gen.GSEA.DF <- function(GSE_file, AnnotFile, GWAS_file){

  #This file returns the the dataframe in the final form for synapse

  #Inputs:
  #GSE_file: List of modules with genes in each
  #AnnotFile: GWAS annotation file
  #GWAS_file: GWAS SNP and p-value file

  #Output:
  #GeneFrame: dataframe containing most enriched module and
  #           pvalue of each gene in the module
  #           also contains p-value of the module computer
  #           by MAGMA

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


Compile.Pval.AllMethods <- function(FolderName){
  #Compile Pvalues of all methods into a common file
  #assumes only relevant output files are in the folder 
  
  FileNames <- list.files(FolderName)
  #print(FileNames)
  
  for (i in 1:length(FileNames)){
    
    tmpStr <- paste(c(FolderName,'/',FileNames[i]),collapse = '')
    temp <- read.table(tmpStr, 
                       skip = 3, header = T)
    
    if (i == 1){
      GWAS_enrich <- temp
    } else {
      GWAS_enrich <- rbind(GWAS_enrich, temp)
    }
    
  }
  
  return(GWAS_enrich)
}


