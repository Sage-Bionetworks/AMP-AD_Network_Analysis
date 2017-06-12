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
Types <- c('consen','megen','metan','speakE','wina')

#5. Visualize modules colored by algorithm type
Net <- GenGraphViz(Dat,pattern, Types)

#6. Visualize modules colored by clusters obtained by 
#graph clustering using markov clustering 
ClustResults <- GenClusteredViz(Net)

#7. Log into Synapse 
synapseClient::synapseLogin()

#Get all the modules in a cluster of interest 
ModClust <- GetModulesInClust(V(Net)$name, 
                              ClustResults$Clusts$Cluster, 1)

#Get all genes in the module set 
ClustGenes <- GetUnionGenes(ModClust)




