setwd('Documents/SageStuff/')

#Defining initial parameters 
FileNames <- list.files('GeneSets/')
Region <- 'DLPFC'
MethodName <- 'consensus'

l <- list()

for (i in 1:length(FileNames)){ 
  
  FileNames[i] <- paste(c('GeneSets/',FileNames[i]),collapse = '')
  }

library(stringr)
source('SumitUtilFcns.R')
source('AMP-AD_Network_Analysis/enrichmentAnalysis/run_amp_ad_enrichment.R')

for (i in 1:length(FileNames)){
  GeneSet <- Create.GeneSet.List(FileNames[i])
  #print('I was here0')
  cat('running gene set',i,'of',length(FileNames),'..\n')
  result <- try({Temp <- run_amp_ad_enrichment_subset(GeneSet, 
                                                      'NA',
                                                      MethodName, 
                                                      Region)
  })

  if (class(result) == "try-error") next; 
  
  #print('I was here')
  
  #tempStr1 <- paste(c('l$mod',i,'$analysis <- Temp'),
  #                 collapse = '')
  tempStr1 <- paste(c('l$\"',
                      str_replace_all(FileNames[i],"[^[:alnum:]]","")
                      , '\" <- Temp'))
  
  #print('I was here2')
  
  #tempStr2 <- paste(c('l$mod',i,'$name <- FileNames[i]'),
  #                  collapse = '')
  
  #print('I was here3')
  
  eval(parse(text = tempStr1))
  #eval(parse(text = tempStr2))
  
  #print('I was here4')
}

saveRDS(l, 'DLPFC_consensus.rds')
