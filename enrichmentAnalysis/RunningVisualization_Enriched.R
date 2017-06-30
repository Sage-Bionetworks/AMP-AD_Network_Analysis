#1. set the working directory and add modules 
setwd('Documents/SageStuff/')
source('GenGraphViz.R')
source('GeneSetAnalysis_Sumit.R')
source('SumitUtilFcns.R')

#2. Load the dataset 
Dat <- read.csv('Job-38889986948603617091033777.csv')
Dat <- data.frame(Dat)

#3. Specify the region of the brain to visualize 
pattern <- "DLPFC"

#4. Specify the algorithms to consider 
Types <- c('consen','megen','metan','speakE','wina')
#Types <- c('consen')


#5. Visualize modules colored by algorithm type
RetGraph <- GenGraphViz(Dat,pattern, Types)
Net <- RetGraph$net 
lt <- RetGraph$l 

#6. Compile the list of all modules at their associated category
FileNames <- list.files('DLPFC_enrichment/')

for (i in 1:length(FileNames)){ 
  
  LoadName <- paste(c('DLPFC_enrichment/',FileNames[i]),collapse='')
  l <- readRDS(LoadName)
  Names <- names(l)
  tmp <- data.frame(AssignEnrichedModule(l,Names[6]))

  if (i == 1){
    Dat <- tmp
  } else {
    Dat <- rbind(Dat, tmp)
  }

}

#Generate color and size profiles for visualization 
LegendNames <- unique(Dat$CatList)
NamesNodes <- names(V(Net))
Labels <- rep(1,length(NamesNodes))
LabelNames <- rep('',length(NamesNodes))
SizeList <- rep(1,length(NamesNodes))

for(i in 1:length(NamesNodes)){
  
  In <- which(Dat$UnqMod %in% NamesNodes[i])
  Labels[i] <- which(LegendNames %in% Dat$CatList[In])
  LabelNames[i] <- Dat$CatList[In]
  SizeList[i] <- -log10(Dat$Pval[In])
  
}

Net2 <- PlotGraphGivenLabels(Net, lt, Labels, SizeList )




