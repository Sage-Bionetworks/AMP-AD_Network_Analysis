synapseClient::synapseLogin()
allMods <- synapseClient::synTableQuery("SELECT * FROM syn9770791")@values
library(dplyr)
fooSummarize <- dplyr::group_by(allMods,brainRegion,ModuleName,method)%>%
  dplyr::summarise(numberOfGenes=length(ModuleName))

fooSummarize2 <- dplyr::group_by(fooSummarize,brainRegion,method) %>%
  dplyr::summarise(numberOfModules = length(method))

barplot(fooSummarize$numberOfGenes,
        log='y',
        col=rainbow(7)[as.factor(fooSummarize$brainRegion)],
        border=NA,
        ylab='Module Size')


png(file='~/Desktop/sizeDistn.png',
    height=800,
    width=1200,
    res=120,
    pointsize = 30)
g <- ggplot2::ggplot(fooSummarize, 
                     ggplot2::aes(x=brainRegion,
                                  y=numberOfGenes,
                                  fill=method))
g <- g + ggplot2::geom_boxplot(position='dodge')
g <- g + ggplot2::scale_y_log10()
g
dev.off()
#splitBr <- lapply(unique(allMods$brainRegion),function(x,y){
#  return(dplyr::filter(y,brainRegion==x))
#},allMods)

#names(splitBr) <- unique(allMods$brainRegion)
#table(splitBr$DLPFC$ModuleName,splitBr$DLPFC$method)



