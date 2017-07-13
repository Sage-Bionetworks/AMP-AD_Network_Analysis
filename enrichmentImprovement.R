synapseClient::synapseLogin()
#consensus

reformatAndCorrectEnrichment <- function(synId){
  afoo <- synapseClient::synGet(synId)
  boo<-readRDS(afoo@filePath) 
  boo2 <- do.call(rbind,boo)
  boo2 <- dplyr::mutate(boo2,adj=p.adjust(boo2$fisherPval,method='fdr'))
  boo2 <- dplyr::filter(boo2,adj<=0.05)
  boo3 <- dplyr::group_by(boo2,ModuleNameFull,geneSet)
  boo4 <- dplyr::summarise(boo3,mean_or=mean(fisherOR),mean_pval=mean(-log10(fisherPval)),nsigcat=length(category))

  return(boo4)
}

reformatAndCorrectEnrichment2 <- function(synId){
  afoo <- synapseClient::synGet(synId)
  boo<-readRDS(afoo@filePath) 
  boo2 <- do.call(rbind,boo)
  boo2 <- dplyr::mutate(boo2,adj=p.adjust(boo2$fisherPval,method='fdr'))
  boo2 <- dplyr::filter(boo2,adj<=0.05)
  boo3 <- dplyr::group_by(boo2,ModuleNameFull,geneSet)
  boo4 <- dplyr::summarise(boo3,mean_or=mean(fisherOR),mean_pval=mean(-log10(fisherPval)),nsigcat=length(category))
  foobar <- list()
  foobar$enrichment <- boo2
  foobar$enrichmentSummary <- boo4
  return(foobar)
}

consensus2 <- reformatAndCorrectEnrichment2('syn10165934')
metanetwork <- reformatAndCorrectEnrichment2('syn10165935')


consensus <- reformatAndCorrectEnrichment('syn10165934')
megena <- reformatAndCorrectEnrichment('syn10166179')
metanetwork <- reformatAndCorrectEnrichment('syn10165935')
rwgcna <- reformatAndCorrectEnrichment('syn10165936')
speakeasy <- reformatAndCorrectEnrichment('syn10165937')
wina <- reformatAndCorrectEnrichment('syn10165938')
summaryDf <- list()
summaryDf$consensus <- consensus
summaryDf$megena <- megena
summaryDf$metanetwork <- metanetwork
summaryDf$rwgcna <- rwgcna
summaryDf$speakeasy <- speakeasy
summaryDf$wina <- wina
addName <- function(x,df){
  df$method <- rep(x,nrow(df))
  return(df)
}
summaryDf2 <- mapply(addName,
                     c('consensus','megena','metanetwork','rwgcna','speakeasy','wina'),
                     summaryDf,
                     SIMPLIFY=FALSE)
summaryDf3 <- do.call(rbind,summaryDf2)

summaryDf3 <- data.frame(summaryDf3,stringsAsFactors = F)
rSynapseUtilities::makeTable(summaryDf3,'enrichment summaries dlpfc July 8 2017','syn2370594')

summaryDf3 <- synapseClient::synTableQuery("SELECT * FROM syn10167932")@values

g <- ggplot2::ggplot(summaryDf3, 
                     ggplot2::aes(x=log10(nsigcat),
                                  color=method,
                                  fill=method))
g <- g + ggplot2::geom_density(alpha = 0.05)
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g

g <- ggplot2::ggplot(summaryDf3, 
                     ggplot2::aes(x=log10(mean_or),
                                  color=method,
                                  fill=method))
g <- g + ggplot2::geom_density(alpha = 0.05)
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g

synapseClient::synapseLogin()
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10163855")@values
library(dplyr)
fooSummarize <- dplyr::group_by(allMods,brainRegion,ModuleName,method)%>%
  dplyr::summarise(numberOfGenes=length(ModuleName))

fooSummarize <- dplyr::mutate(fooSummarize,ModuleNameFull = paste0(ModuleName,brainRegion))

summaryDf3 <- dplyr::left_join(summaryDf3,fooSummarize,by='ModuleNameFull')
summaryDf3$mean_or[!is.finite(summaryDf3$mean_or)] <- NA
summaryDf4<-na.omit(summaryDf3)
summaryDf4 <- dplyr::mutate(summaryDf4,meanORresid = MASS::rlm(log(mean_or) ~ log(numberOfGenes))$resid)
summaryDf4 <- dplyr::mutate(summaryDf4,meanORresidexp=exp(meanORresid))
summaryDf4 <- dplyr::mutate(summaryDf4,ncatresid = MASS::rlm(log(nsigcat) ~ log(numberOfGenes))$resid)


g <- ggplot2::ggplot(summaryDf4, 
                     ggplot2::aes(x=method.x,
                                  y=ncatresid,
                                  color=method.x,
                                  fill=method.x))
g <- g + ggplot2::geom_boxplot(alpha = 0.2)
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g

g <- ggplot2::ggplot(summaryDf4, 
                     ggplot2::aes(x=method.x,
                                  y=log(mean_pval),
                                  color=method.x,
                                  fill=method.x))
g <- g + ggplot2::geom_boxplot(alpha = 0.2)
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g


g <- ggplot2::ggplot(summaryDf4, 
                     ggplot2::aes(x=method.x,
                                  y=meanORresid,
                                  color=method.x,
                                  fill=method.x))
g <- g + ggplot2::geom_boxplot(alpha = 0.2)
#g <- g + ggplot2::scale_y_log10()
g <- g + ggplot2::theme_grey(base_size = 20) 
g

means1 <- dplyr::group_by(summaryDf4,method.x) %>% dplyr::summarize(meanncat=median(ncatresid),meanor=median(meanORresid))
