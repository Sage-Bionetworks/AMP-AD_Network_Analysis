##load top targets
topTargets <- readLines('amp_ad_top_targets')

##load network results
load('aggregate_module_mainfest.rda')


synapseClient::synapseLogin()

##load deg lists
degResObj <- synapseClient::synGet("syn10496554")
load(degResObj@filePath)

degDf <- utilityFunctions::list2df(amp.ad.de.geneSets)
degDfTop <- dplyr::filter(degDf,value%in%topTargets)
wup<-grep('UP',degDfTop$key)
degDfTop$direction <- rep('DOWN',nrow(degDfTop))
degDfTop$direction[wup] <- 'UP'
View(degDfTop)
getBr <- function(x){
  return(strsplit(x,'\\.')[[1]][1])
}
degDfTop$br <- sapply(degDfTop$key,getBr)
#degDfTop <- dplyr::mutate(degDfTop,br = getBr(key))
g <- ggplot2::ggplot(degDfTop, 
                     ggplot2::aes(value,
                                  key,
                                  shape = direction,
                                  color = br))
g <- g + ggplot2::geom_count()
g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
g

###summarize things

fob <- dplyr::group_by(degDfTop,value,br)
bar <- dplyr::summarise(fob,count = length(unique(br)))

bar2 <- dplyr::group_by(bar,value)
bar3 <- dplyr::summarise(bar2,count = sum(value!=0))
#bar3 <- dplyr::arrange(bar3,desc(count))
bar3$value <- factor(bar3$value,levels=bar3$value)
g <- ggplot2::ggplot(bar3,
                     ggplot2::aes(value,
                                  count))
g <- g + ggplot2::geom_col()
g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
g <- g + ggplot2::xlab('Top Targets')
g <- g + ggplot2::ylab('# Brain Regions differentially expressed')
g


####get module assignments

constructQuery <- function(x,synId){
  str <- paste0('SELECT * FROM ',synId,' WHERE ')
  str2 <- paste0(x,collapse = '\' OR external_gene_name = \'')
  str2 <- paste0(str,'external_gene_name = \'',str2)
  str2 <- paste0(str2,'\'')
  return(str2)
}

queryStr <- constructQuery(topTargets,'syn11182793')
modTable <- synapseClient::synTableQuery(queryStr)@values



#get hub score for each module x gene

get_hub_score <- function(gene,module,masterList){
  return(dplyr::filter(masterList[[module]]$anno,external_gene_name==gene)$hubs)
}
modTable$hubScore <- mapply(get_hub_score,
                            modTable$external_gene_name,
                            modTable$ModuleNameFull,
                            MoreArgs = list(masterList=fullManifest))

g <- ggplot2::ggplot(modTable, 
                     ggplot2::aes(external_gene_name,
                                  ModuleName,
                                  size=hubScore))
g <- g + ggplot2::geom_count()
g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
g

bar5 <- dplyr::group_by(modTable,external_gene_name)
bar6 <- dplyr::summarise(bar5,medianScore=median(hubScore))
bar6 <- dplyr::arrange(bar6,desc(medianScore))
bar6$external_gene_name <- factor(bar6$external_gene_name,levels=bar6$external_gene_name)
g <- ggplot2::ggplot(bar6,
                     ggplot2::aes(external_gene_name,
                                  medianScore))
g <- g + ggplot2::geom_col()
g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
g <- g + ggplot2::xlab('Top Targets')
g <- g + ggplot2::ylab('Median module connectivity')
g

foobar <- dplyr::left_join(bar6,
                           bar3,
                           by=c('external_gene_name'='value'))

####cross consortia top targets based on same metrics
###get network degrees

getHubs <- function(modbr){
  bicNets <- synapseClient::synTableQuery(paste0("SELECT * FROM syn8681664 where ( (method = \'bic\') and (tissueTypeAbrv = \'",modbr,"\' )  and ( assay = \'RNAseq\'))"))@values
  load(synapseClient::synGet(bicNets$id[1])@filePath)
  library(Matrix)
  res <- list()
  res$adjacencyMatrix <- as.matrix(bicNetworks$network)
  hubs <- data.frame(hubs = rowSums(res$adjacencyMatrix+t(res$adjacencyMatrix)),GeneID=rownames(res$adjacencyMatrix),brainRegion = rep(modbr,nrow(res$adjacencyMatrix)),stringsAsFactors=F)
  
  
  return(hubs)
}

hubDf <- lapply(c('DLPFC','STG','IFG','FP','PHG','TCX','CBE'),getHubs)
hubDf <- do.call('rbind',hubDf)

map1 <- utilityFunctions::convertEnsemblToHgnc(unique(hubDf$GeneID))
hubDf <- dplyr::left_join(hubDf,map1,by=c('GeneID'='ensembl_gene_id'))
