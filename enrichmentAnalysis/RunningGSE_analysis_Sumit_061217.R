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
pattern <- "CER"

#4. Specify the algorithms to consider 
Types <- c('consen','megen','metan','speakE','wina')

#5. Visualize modules colored by algorithm type
RetGraph <- GenGraphViz(Dat,pattern, Types)
Net <- RetGraph$net 
l <- RetGraph$l 

#6. Visualize modules colored by clusters obtained by 
#graph clustering using markov clustering 
ClustResults <- GenClusteredViz(Net,l)

#7. Log into Synapse 
synapseClient::synapseLogin()

#8. Get all the modules in a cluster of interest 
ModClust <- GetModulesInClust(V(Net)$name, 
                              ClustResults$Clusts$Cluster, 1)

#9. Get all genes in the module set 
ClustGenes <- GetUnionGenes(ModClust)

#10. Get all module names for a particular module type 
PatternMod <- 'consen'
In1 <- grep(PatternMod,V(Net)$name)
ModNames <- V(Net)$name[In1]

#11. CreateFile Gene set file for MAGMA 
source('CreateMagmaFiles.R')
OutputFileName <- 'MagmaModuleFile.txt'
Create.MAGMA.GeneLists(ModNames, 
                       OutputFileName = OutputFileName)



