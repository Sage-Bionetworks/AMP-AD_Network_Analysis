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
View(rosmapClinical)


####read in uncensored age file
rosmapUncensoredAges <- data.table::fread(synapseClient::getFileLocation(rosmapUncensoredAgesObj),
                                          data.table=F)
View(rosmapUncensoredAges)

####merge uncensored ages with clinical file
rosmapClinical <- dplyr::left_join(rosmapClinical,rosmapUncensoredAges,by='projid')
View(rosmapClinical)

###cognitive decline
rosmapCogDec1 <- data.table::fread(synapseClient::getFileLocation(rosmapCogDecline1Obj),
                                   data.table=F)
rosmapCogDec1 <- dplyr::select(rosmapCogDec1,-Sample)
View(rosmapCogDec1)
rosmapClinical <- dplyr::left_join(rosmapClinical,
                                   rosmapCogDec1,
                                   by=c("projid"="ProjectID"))

rosmapClinical$apoe_genotype <- as.factor(rosmapClinical$apoe_genotype)
rosmapClinical$ceradsc <- as.factor(rosmapClinical$ceradsc)
rosmapClinical$braaksc <- as.factor(rosmapClinical$braaksc)
rosmapClinical$cogdx <- as.factor(rosmapClinical$cogdx)
####simple tests
lmObj <- lm(cogn_global_slope ~ age_death+ceradsc+braaksc+cogdx+apoe_genotype,data = rosmapClinical)
summary(lmObj)

library(dplyr)
KEY <- read.csv(synapseClient::getFileLocation(rosmapIdMapObj))%>% 
  dplyr::filter(mrna_data == 1) %>%
  dplyr::select(projid, mrna_id) %>%
  tidyr::separate(mrna_id, c('a','b','batch'), sep = '_') %>% 
  tidyr::unite(Sampleid, a, b) %>%
  dplyr::select(-batch) %>%
  unique

expressionDataObj <- synapseClient::synGet('syn8456638')
