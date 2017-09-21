#Log into Synapse 
synapseClient::synapseLogin()

#Select brain region
brainRegion <- 'STG'

#Select manifestId where all module information is stored 
manifestId <- 'syn10309369'
seed <- 5

#Getting the module information 
queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",brainRegion,"' ) )")
indMods <- synapseClient::synTableQuery(queryString)@values

#Converting modules to right format for consensus 
source('AMP-AD_Network_Analysis/UtilityFiles/ConvertsModsToMat.R')
InstanceList <- ConvertModsToMat(indMods)

#loading packages required for consensus module generation 
library(AggMethods)
setwd('AggMethods/')
library('devtools')
library('roxygen2')
document()
update.packages('AggMethods')
library(AggMethods)

#Perform consensus module generation using CC-Pivot
G <- Corr2Cons_Mod(InstanceList)

Ep <- G$Ep
Em <- G$Em
n<- sqrt(length(Ep))

Lab <- rep(0,n)
document()
update.packages('AggMethods')
CC.Pivot.Global(G$V,1)

#If required perform the following step to generate module using CombinedClusteringWithReps
##Cons <- CombinedClusteringWithRepsMod(InstanceList, 25)