#####grab degs
synapseClient::synapseLogin()
degResObj <- synapseClient::synGet("syn10163525")
degRes <- data.table::fread(degResObj@filePath,
                            data.table=FALSE)
degRes <- dplyr::filter(degRes,
                        adj.P.Val<=0.05)
degRes <- dplyr::mutate(degRes,Comparison2 = paste0(Comparison,'_',Direction))

source('enrichmentAnalysis/run_amp_ad_enrichment.R')

listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}

degList <- lapply(unique(degRes$Comparison2),
                           listify,
                           degRes$Gene.ID,
                           degRes$Comparison2)
names(degList) <- unique(degRes$Comparison2)
reformatNames <- function(x){
  foo1<-strsplit(x,'\\.')
  foo2<-sapply(foo1,function(y) y[1])
  return(foo2)
}
degList <- lapply(degList,reformatNames)

degResults <- run_amp_ad_enrichment(degList,
                                    "degs",
                                    hgnc=FALSE)

###reorganize deg results
#split off brain region
parseDegName <- function(x){
  library(dplyr)
  foo1<-strsplit(x,'_')[[1]]
  br<-foo1[1]
  dir <- foo1[length(foo1)]
  cate <- paste0(foo1[2:(length(foo1)-1)],collapse='_')
  if(length(grep(paste0('.',br,'_'),cate))>0){
    cate<-gsub(paste0('.',br,'_'),'.',cate)
  }
  
  c('brainRegion'=br,
    'Direction'=dir,
    'reducedCategory'=cate,
    'Category'=x) %>% return
}
categoryKey <- sapply(unique(degResults$category),
                      parseDegName)
categoryKey <- t(categoryKey)
categoryKey <- data.frame(categoryKey,stringsAsFactors=F)

degResults2 <- dplyr::left_join(degResults,categoryKey,by=c('category'='Category'))

degResults2 <- dplyr::mutate(degResults2,Z=qnorm(fisherPval,lower.tail=F))



summaryDegManifest <- dplyr::group_by(degResults2,ModuleNameFull,Direction,reducedCategory) %>% 
  dplyr::summarise(medianZ = median(Z))

g <- ggplot2::ggplot(summaryDegManifest, 
                     ggplot2::aes(x=Direction,
                                  y=medianZ,
                                  fill=reducedCategory))
g <- g + ggplot2::geom_boxplot(position='dodge')
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g

#dplyr::summarise(numberOfGenes=length(ModuleName)

#categoryKey <- categoryKey[!duplicated(categoryKey),]
#####MAGMA results
magmaResults <- synapseClient::synTableQuery("SELECT * FROM syn10380432")@values

#####cell type results
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
cellMarkers <- GeneSets$Cell_Markers
cellTypeResults <- run_amp_ad_enrichment(cellMarkers,
                                        "celltypes",
                                        hgnc=TRUE)
#####combined manifest
fullManifest <- rbind(degResults,
                      cellTypeResults)

magmaReformat <- dplyr::select(magmaResults,SET,P)

colnames(magmaReformat) <- c('ModuleNameFull','magmaPval')
magmaReformat <- dplyr::mutate(magmaReformat,magmaZ=qnorm(magmaPval,lower.tail=F))

summaryDegManifest2 <- dplyr::left_join(summaryDegManifest,magmaReformat)
summaryDegManifest2 <- dplyr::mutate(summaryDegManifest2,combZ=magmaZ/2+medianZ/2)
summaryDegManifest2 <- dplyr::arrange(summaryDegManifest2,desc(medianZ))
#split by each category
fxn1 <- function(x,y){
  foobar <- dplyr::filter(y,reducedCategory==x)
  return(foobar)
}
splitSummaries <- lapply(unique(summaryDegManifest2$reducedCategory),fxn1,summaryDegManifest2)
names(splitSummaries) <- unique(summaryDegManifest2$reducedCategory)

##just take top from each up/down
fxn2 <- function(x){
  foobar1 <- dplyr::filter(x,Direction=='DOWN')
  foobar2 <- dplyr::filter(x,Direction=='UP')
  return(c('down_mod'=foobar1$ModuleNameFull[1],
           'up_mod'=foobar2$ModuleNameFull[1]))
}
getMods <- sapply(splitSummaries,
                  fxn2)
topMods <- t(getMods)

#magmaReformat$category <- rep('MAGMA',nrow(magmaReformat))
#magmaReformat$fisherOR <- rep(NA,nrow(magmaReformat))
fullManifest <- dplyr::select(fullManifest,ModuleNameFull,category,fisherPval,fisherOR)
fullManifest <- rbind(fullManifest,magmaReformat)
fullManifest <- dplyr::mutate(fullManifest,Z = qnorm(fisherPval,lower.tail=F))
fullManifestSquare <- dplyr::select(fullManifest,ModuleNameFull,category,Z)
fullManifestSquare <- tidyr::spread(fullManifestSquare,ModuleNameFull,Z)
rownames(fullManifestSquare) <- fullManifestSquare$category
fullManifestSquare <- dplyr::select(fullManifestSquare, -category)
fullManifestSquare <- data.matrix(fullManifestSquare)
fullManifestSquare[!is.finite(fullManifestSquare)] <- NA
fullManifestSquare <- t(fullManifestSquare)
#fullManifestSquare[is.na(fullManifestSquare)] <- 0
foobar <- apply(fullManifestSquare,1,median,na.rm=T)
#####combined score


#####top modules
