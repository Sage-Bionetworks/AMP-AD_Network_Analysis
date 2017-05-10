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



g <- ggplot2::ggplot(fooSummarize, ggplot2::aes(numberOfGenes,fill=method))
g <- g + ggplot2::geom_bar(position="dodge")
g
#splitBr <- lapply(unique(allMods$brainRegion),function(x,y){
#  return(dplyr::filter(y,brainRegion==x))
#},allMods)

#names(splitBr) <- unique(allMods$brainRegion)
#table(splitBr$DLPFC$ModuleName,splitBr$DLPFC$method)



