
load('aggregate_module_mainfest.rda')
gwasVariants <- data.table::fread("ROSMAP_Mayo_Extracted_SNPs_forBen.txt",data.table=F)

getTop5pcs <- function(x){
  pcmat <- svd(scale(x$expr))$u[,1:5]
  rownames(pcmat) <- rownames(x$expr)
  return(pcmat)
}
pcmats <- lapply(fullManifest,getTop5pcs)

###get rosmap individual mapping function
rosmapIdMapObj <- synapseClient::synGet('syn3382527')
library(dplyr)
rosmapIDmap <- data.table::fread(rosmapIdMapObj@filePath,data.table=F)




KEY <- read.csv(synapseClient::getFileLocation(rosmapIdMapObj))%>%
  dplyr::filter(mrna_data == 1) %>%
  dplyr::select(projid, mrna_id) %>%
  tidyr::separate(mrna_id, c('a','b','batch'), sep = '_') %>%
  tidyr::unite(Sampleid, a, b) %>%
  dplyr::select(-batch) %>%
  unique

View(KEY)

KEY2 <- dplyr::left_join(KEY,dplyr::select(rosmapIDmap,projid,gwas_id))

###test it out on one rosmap module
rosmapMod <- pcmats[[13]]

###pull eigengenes from all modules


colnames(rosmapMod) <- paste0('pc',1:5)
rosmapMod <- data.frame(rosmapMod,stringsAsFactors = F)
rosmapMod$Sampleid <- rownames(rosmapMod)
rosmapMod <- dplyr::left_join(rosmapMod,KEY2)
rosmapMod <- dplyr::left_join(rosmapMod,gwasVariants,by=c('gwas_id'='ID'))

rosmapMod2 <- na.omit(rosmapMod)

summary(lm(pc1 ~ rs3865444,data=rosmapMod))

lm_wrapper <- function(y,x){
  return(summary(lm(y~x))$coef[2,4])
}

gwas_test <- utilityFunctions::outerApplyParallel(lm_wrapper,
                                     dplyr::select(rosmapMod,dplyr::starts_with('rs')),
                                     dplyr::select(rosmapMod,dplyr::starts_with('pc')))

gap::qqunif(c(gwas_test))
pheatmap::pheatmap(-log10(gwas_test))

###get mayo individual mapping function
View(fullManifest[[33]])
#fullManifest[[33]]$expr[1:5,1:5]
mayoMod <- pcmats[[5]]
colnames(mayoMod) <- paste0('pc',1:5)
mayoMod <- data.frame(mayoMod,stringsAsFactors=F)
mayoMod$ID <- unlist(strsplit(rownames(mayoMod),'_CER'))
mayoMod <- dplyr::left_join(mayoMod,gwasVariants)
gwas_test <- utilityFunctions::outerApplyParallel(lm_wrapper,
                                                  dplyr::select(mayoMod,dplyr::starts_with('rs')),
                                                  dplyr::select(mayoMod,dplyr::starts_with('pc')))

gap::qqunif(c(gwas_test))
pheatmap::pheatmap(-log10(gwas_test))



#####
load(synapseClient::synGet('syn10737440')@filePath)

rosmapEG <- fullEigengeneSet$rosmapDLPFC
rosmapEG<-lapply(rosmapEG,function(x,rn){
  if(!is.na(x)){
    rownames(x) <- rn
    colnames(x) <- paste0('pc',1:5)
  }
  return(x)
},rownames(pcmats[[13]]))
#dim(pcmats[[13]])

lm_wrapper <- function(y,x){
  return(summary(lm(y~x))$coef[2,4])
}

run_module_gwas_analysis <- function(modPC,npc,KEY2,gwasVariants){
  if(!is.na(modPC)){
    rosmapMod <- data.frame(modPC,stringsAsFactors = F)
    rosmapMod$Sampleid <- rownames(rosmapMod)
    rosmapMod <- dplyr::left_join(rosmapMod,KEY2)
    rosmapMod <- dplyr::left_join(rosmapMod,gwasVariants,by=c('gwas_id'='ID'))
    rosmapMod2 <- na.omit(rosmapMod)
    gwas_test <- utilityFunctions::outerApply(lm_wrapper,
                                                      dplyr::select(rosmapMod2,dplyr::starts_with('rs')),
                                                      dplyr::select(rosmapMod2,dplyr::starts_with('pc')))
    if(min(gwas_test) < 0.05/115){
      print(npc)
      print(min(gwas_test))
      a1 <- which(gwas_test == min(gwas_test),T)
      print(rownames(gwas_test)[a1[1]])
      print(colnames(gwas_test)[a1[2]])
      #return(gwas_test)
    }
    return(gwas_test)
  }else{
    return(NA)
  }
  #return(gwas_test)
  
}

rosmapEG2 <- rosmapEG[grep('DLPFC',names(rosmapEG))]
aaa<-mapply(run_module_gwas_analysis,
       rosmapEG2,
       names(rosmapEG2),
       MoreArgs = list(KEY2=KEY2,
                       gwasVariants=gwasVariants),
       SIMPLIFY = F)


pheatmap::pheatmap(-log10(aaa$megenac1_67DLPFC))


rosmapMod <- data.frame(rosmapEG2$megenac1_67DLPFC,stringsAsFactors = F)
rosmapMod$Sampleid <- rownames(rosmapMod)
rosmapMod <- dplyr::left_join(rosmapMod,KEY2)
rosmapMod <- dplyr::left_join(rosmapMod,gwasVariants,by=c('gwas_id'='ID'))
rosmapMod2 <- na.omit(rosmapMod)


###
synapseClient::synapseLogin()
rosmapClinicalObj <- synapseClient::synGet('syn3191087')
rosmapUncensoredAgesObj <- synapseClient::synGet('syn7116000')
rosmapIdMapObj <- synapseClient::synGet('syn3382527')
rosmapCogDecline1Obj <- synapseClient::synGet('syn6182375')
rosmapCogDecline2Obj <- synapseClient::synGet('syn6182376')
rosmapMetanetworkObj <- synapseClient::synGet('syn8268669')

####thanneer's code for fixing ids
#https://github.com/th1vairam/Brain_Reg_Net/blob/ad12b544d6d6c23be2bf94eed53a6b4f75d154d1/code/Rmd/ROSMAP_REPROCESSED.Rmd

####read in clinical file
rosmapClinical <- data.table::fread(synapseClient::getFileLocation(rosmapClinicalObj),
                                    data.table=F)
rosmapClinical <- dplyr::select(rosmapClinical,-V1,-age_death,-age_first_ad_dx,-age_at_visit_max)
#View(rosmapClinical)


####read in uncensored age file
rosmapUncensoredAges <- data.table::fread(synapseClient::getFileLocation(rosmapUncensoredAgesObj),
                                          data.table=F)
#View(rosmapUncensoredAges)

####merge uncensored ages with clinical file
rosmapClinical <- dplyr::left_join(rosmapClinical,rosmapUncensoredAges,by='projid')
#View(rosmapClinical)

###cognitive decline
rosmapCogDec1 <- data.table::fread(synapseClient::getFileLocation(rosmapCogDecline1Obj),
                                   data.table=F)
rosmapCogDec1 <- dplyr::select(rosmapCogDec1,-Sample)

rosmapCogDec2 <- data.table::fread(synapseClient::getFileLocation(rosmapCogDecline2Obj),
                                   data.table=F)

mungeIds <- function(x){
  #ros or map
  isROS<-grep('ROS',x)
  isMAP<-grep('MAP',x)
  if(length(isROS)>0){
    rosId<-strsplit(x,'ROS')[[1]][2]
    return(as.numeric(rosId))
  }else if(length(isMAP)>0){
    mapId <- strsplit(x,'MAP')[[1]][2]
    return(as.numeric(mapId))
  }else{
    return(as.numeric(x))
  }
}
projids2 <- sapply(rosmapCogDec2$CollaboratorParticipantId,mungeIds)
rosmapCogDec2$ProjectID <- projids2
rosmapCogDec2 <- dplyr::select(rosmapCogDec2,ProjectID,cogn_global_slope)
rosmapCogDec <- rbind(rosmapCogDec1,rosmapCogDec2)

#View(rosmapCogDec1)
rosmapClinical <- dplyr::left_join(rosmapClinical,
                                   rosmapCogDec,
                                   by=c("projid"="ProjectID"))
rosmapClinical$apoe_genotype <- as.factor(rosmapClinical$apoe_genotype)
rosmapClinical$ceradsc <- as.factor(rosmapClinical$ceradsc)
rosmapClinical$braaksc <- as.factor(rosmapClinical$braaksc)
rosmapClinical$cogdx <- as.factor(rosmapClinical$cogdx)

rosmapMod2 <- dplyr::left_join(rosmapMod2,
                               rosmapClinical,by='projid')

gap::qqunif(c(gwas_test))