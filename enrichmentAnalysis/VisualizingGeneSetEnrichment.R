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

#7. Compile all module results for a enrichment category 
GeneSetName <- 'BioCarta_2013'
for (i in 1:length(Types)){
  tmp <- data.frame(EnrList[[Types[i]]][[GeneSetName]])
  
  if (i == 1){
    ResFrame <- tmp
  } else {
    ResFrame <- rbind(ResFrame, tmp)
  }
}

source('SumitUtilFcns.R')
source('CreateMagmaFiles.R')

#8. For each module obtain the lowest p-value and most enriched 
#   gene-set 
EnrByCat <- AssignEnrichedModule(ResFrame)


#9. For each module in the graph, obtain the MAGMA p-value
Mod_pval <- rep(1,length(names(V(Net))))
MAGMA_pval <- 
  readRDS('GWAS_files/DLPFC_results_magma/Magma_AllMethods.Rda')
In1 <- which(names(V(Net)) %in% MAGMA_pval$SET)
Mod_pval[In1] <- MAGMA_pval$P


#10. For each module obtain the geneset type 
LegendNames <- unique(EnrByCat$CatList)
Labels <- rep(1,length(Mod_pval))
for (i in 1:length(Mod_pval)){
  Labels[i] <- which( LegendNames %in% EnrByCat$CatList[i])
}


#11. Create visualization 
Net2 <- PlotGraphGivenLabels(Net, l, Labels, -log10(Mod_pval),
                             sizeMin = .5)


library(ggplot2)
#12. Create barplot for a particular module 
ModType <- 'consensus'
ModNameEnr <- 'consensus37DLPFC'
TempEnr <- EnrList[[ModType]][[GeneSetName]]
In <- which(TempEnr$ModuleNameFull %in% ModNameEnr)
tmp <- list()
Category<- TempEnr$category[In]
Log10_pval <- -log10(TempEnr$fisherPval[In])
I = sort(Log10_pval, decreasing = T, index.return = T)
I <- I$ix
tmp$Category<- Category[I]
tmp$Log10_pval<- Log10_pval[I]
tmp <- data.frame(tmp, stringsAsFactors = F)
p <- ggplot(data = tmp, aes(x = Category, y = Log10_pval)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, size = 5)) + 
  scale_x_discrete(limits = tmp$Category)







