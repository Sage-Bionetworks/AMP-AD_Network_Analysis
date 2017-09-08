#This script explains how to run the GenGraphViz() and 
#GeneSetAnalysis_Sumit() function modules 

#1. set the working directory and add modules 
setwd('Documents/SageStuff/')
source('GenGraphViz.R')
source('GeneSetAnalysis_Sumit.R')

#2. Load the module dataset dataset 
Dat <- read.csv('Job-38889986948603617091033777.csv')
Dat <- data.frame(Dat)

Dat$fisherOR <- log(1 + Dat$fisherOR)

#3. Specify the region of the brain to visualize 
pattern <- "DLPFC"

#4. Specify the algorithms to consider 
Types <- c('consen','megen','metan','speakE','wina')
#Types <- c('consen')


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
PatternMod <- 'megen'
In1 <- grep(PatternMod,names(V(Net)))
ModNames <- V(Net)$name[In1]

#11. CreateFile Gene set file for MAGMA 
source('CreateMagmaFiles.R')
OutputFileName <- 
  './GWAS_files/DLPFC_genesets_0705/MagmaModuleFile_Consensus.txt'
Create.MAGMA.GeneLists(ModNames, 
                       OutputFileName = OutputFileName, 
                       manifestId = 'syn10146524')

#12. Plot histogram of pValues 
GWAS_enrich <- read.table('ModuleGSEA_DLPFC.log.sets.out', 
                          skip = 3, header = T)
xlab <- 'Modules'
ylab <- '-log10(pval)'
barplot(-log10(GWAS_enrich$P), xlab = xlab, ylab = ylab)

#13. Find the most enriched module 
Min_In <- which.min(GWAS_enrich$P)
Min_mod <- as.vector(droplevels(GWAS_enrich$FULL_NAME[Min_In]))

#14. Get minimum SNP p-values for genes in this module 
GenePvalList <- GWAS.Enrich.Modules(Min_mod, 
                                    AnnotFile = 'testOut.genes.annot'
                                    , GWAS_file = 'IGAP_stage_1.txt',
                                    manifestId = 'syn10158502')

#15. Find significant genes in the module 
In_imp <- which(GenePvalList$Pval_min < 1e-5)
GenePvalList$Genes[In_imp]

#16. plotting bargraph of pvalues 
barplot(-log10(GenePvalList$Pval_min))

#17. Write to Rda file 
GenePvalList$MagmaPval <- 
  rep(GWAS_enrich$P[Min_In],length(GenePvalList$Pval_min))
GenePvalList$ModuleName <- 
  rep(Min_mod,length(GenePvalList$Pval_min))
GeneFrame <- data.frame(GenePvalList)
OutFile <- paste(c('./GeneSetResults/',Min_mod,'_GeneList.Rda'),collapse = '')
save(GeneFrame, file = OutFile)

  
  

