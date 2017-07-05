#load the file containing all the pathway definitions 

synapseClient::synapseLogin()
synId <- 'syn4867851'

foo <- synapseClient::synGet(synId)
fileLoc <- synapseClient::getFileLocation(foo)
load(fileLoc)
SetNames <- names(GeneSets)

#Run enrichment analysis 

# l <- list()
# for (i in 1:length(SetNames)){
#   tmp <- paste(c('l$\'',SetNames[i],'\' <- c(1,2,3)'), collapse = '')
#   #l$SetNames[i] <- c(1,2,3)
#   eval(parse(text = tmp))
# }

#Run enrichment analysis 
#Types <- c('consensus','megena','metanetwork','speakEasy','wina')
Types <- c('consensus')
Region <- 'DLPFC'

manifestId <- 'syn10158502'

source('AMP-AD_Network_Analysis/enrichmentAnalysis/run_amp_ad_enrichment.R')

l <- list()

for (i in 1:length(Types)){
  
  cat('Module type = ',Types[i],' ...\n')
  
  l2 <- list()
  
  j <- 1
  ctr <- 1
  ctr_max <- 5
  
  while (j <= length(SetNames)){
    
    cat('Running enrichment for',SetNames[j], '...\n')
    
    Temp <- c()
    result <- try({Temp <- run_amp_ad_enrichment_subset(GeneSets[[SetNames[j]]],
                                                        geneSetName = SetNames[j],
                                                        method = Types[i],
                                                        brainRegion = Region,
                                                        manifestId = manifestId)})
    if((class(result) == "try-error") & length(Temp)==0 ){
      if (ctr <= ctr_max){
        cat('Failed trial no.',ctr,'...\n')
        ctr <- ctr + 1
        next
      } else {
        ctr <- 1
        cat('Giving up ... \n')
        j <- j + 1 
        next 
      } 
      
    }
    
    
    tmp <- paste(c('l2$\'',SetNames[j],'\' <- Temp'), collapse = '')
    eval(parse(text = tmp))
    #head(l2)
    j <- j + 1
  }
  
  tmp2 <- paste(c('l$\'',Types[i],'\' <- l2'), collapse = '')
  eval(parse(text = tmp2))
  head(l)
  
}


saveRDS(l, "DLPFC_consen_0705.rds")  



