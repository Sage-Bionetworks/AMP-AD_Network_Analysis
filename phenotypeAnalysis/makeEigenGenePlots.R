synapseClient::synapseLogin()

moduleManifestCombined4 <- synapseClient::synTableQuery("SELECT * FROM syn9770842")@values


#####stratify p-values by 1) brain Region of association, 2) brain region of module, 3) method

#eigengene1
png(file='~/Desktop/eigenGenePlot.png',
    height=800,
    width=1200,
    res=120,
    pointsize = 30)
g <- ggplot2::ggplot(moduleManifestCombined4, 
                     ggplot2::aes(x=brainRegion,
                                  y=-log10(eigengeneAggregate),
                                  fill=brainRegionAssociation))
g <- g + ggplot2::geom_boxplot(position='dodge')
g <- g + ggplot2::labs( x = 'Brain Region Module Inferred in',
               y = '-log10 p-value',
               title= 'Association of module expression with AD/Control Status')
g
dev.off()
