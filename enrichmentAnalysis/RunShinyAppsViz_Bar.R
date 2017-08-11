#This script explains how to run the GenGraphViz() and 
#GeneSetAnalysis_Sumit() function modules 

#1. set the working directory and add modules 
setwd('Documents/SageStuff/')
source('GenGraphViz.R')
source('GeneSetAnalysis_Sumit.R')

#2. Load the dataset 
Dat <- read.csv('Job-38889986948603617091033777.csv')
Dat <- data.frame(Dat)

#3. Specify the region of the brain to visualize 
pattern <- "DLPFC"

#4. Specify the algorithms to consider 
Types <- c('consensus','megena','metanetwork','speakEasy','wina')
#Types <- c('consen')


#5. Visualize modules colored by algorithm type
RetGraph <- GenGraphViz(Dat,pattern, Types)
Net <- RetGraph$net 
l <- RetGraph$l 

#6. Read the list containing enrichment for all modules 
EnrList <- readRDS('DLPFC_allMods_0626.rds')

source('SumitUtilFcns.R')
source('CreateMagmaFiles.R')
library(shiny)

#run network visualization Shiny app 
runApp('ShinyNetViz/')

#run barplot visualization Shiny app
runApp('ShinyEnrichBar/')







