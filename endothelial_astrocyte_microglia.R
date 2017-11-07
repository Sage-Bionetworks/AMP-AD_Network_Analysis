load('aggregate_module_mainfest.rda')

#CBE Turquoise, PHG Turquoise, STG Turquoise, TCX turquoise, IFG Turquoise, DLPFC blue, TCX blue

mods <- c("aggregateCBEturquoiseCBE",
          "aggregatePHGturquoisePHG",
          "aggregateSTGturquoiseSTG",
          "aggregateTCXturquoiseTCX",
          "aggregateIFGturquoiseIFG",
          "aggregateDLPFCblueDLPFC",
          "aggregateTCXblueTCX")
manifest <- fullManifest[mods]

edgeLists <- lapply(manifest,function(x){
  library(Matrix)
  foobar <- utilityFunctions::convertAdjacencyToEdgeList(x$adjacencyMatrix)
  foobar1 <- utilityFunctions::convertEnsemblToHgnc(unique(c(foobar)))
  foobar <- data.frame(foobar,stringsAsFactors=F)
  colnames(foobar) <- c('node1','node2')
  foobar <- dplyr::left_join(foobar,foobar1,by=c('node1'='ensembl_gene_id'))
  foobar <- dplyr::left_join(foobar,foobar1,by=c('node2'='ensembl_gene_id'))
  return(foobar)
})

filenames <- c('cbe_turquoise.csv',
               'phg_turquoise.csv',
               'stg_turquoise.csv',
               'tcx_turquoise.csv',
               'ifg_turquoise.csv',
               'dlpfc_turquoise.csv',
               'tcx_blue.csv')
mapply(function(x,y){
  write.csv(x[,3:4],file=y,quote=F)
},edgeLists,filenames)

#make table of cell types
synapseClient::synapseLogin()
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
cellTypeDf <- utilityFunctions::list2df(GeneSets$Cell_Markers)
write.csv(cellTypeDf,file='cellTypeAnnotations.csv',quote=F)
#GeneSets$Cell_Markers


#work on expression plots

#plot that shows summary of how much genes are up or down across the 7 modules for each deg set

#bar plot for each deg fold change estimate sorted by effect size colored by direction and significance for each pair of module x deg_set

#get fold changes from synapse
library(dplyr)
synapseClient::synapseLogin()
degResProv<-synapseClient::synGetActivity('syn10496554')
synapseIdsOfDegs <- sapply(degResProv@properties$used,function(x) x$reference$targetId) %>% unlist
sapply(synapseIdsOfDegs,synapseClient::onWeb)

#load deg results into a list of tables
degResultTable <- rSynapseUtilities::loadDelimIntoList(synapseIdsOfDegs)

###DO EACH BRAIN REGION then concatenate

#model,gene,logFC,sig,brainRegion
DLPFCsummary <- dplyr::select(degResultTable$syn8456721,
                              Model,
                              Comparison,
                              logFC,
                              adj.P.Val,
                              ensembl_gene_id,
                              hgnc_symbol,
                              Region)

DLPFCsummary$model <- paste0(DLPFCsummary$Model,'.',DLPFCsummary$Comparison)
DLPFCsummary <- dplyr::select(DLPFCsummary,-Model,-Comparison)
DLPFCsummary <- dplyr::mutate(DLPFCsummary,significant = adj.P.Val <= 0.05)
DLPFCsummary <- dplyr::mutate(DLPFCsummary,gene = ensembl_gene_id)
DLPFCsummary <- dplyr::select(DLPFCsummary,-adj.P.Val)

#split into list for each model
DLPFCsummaryList <- plyr::dlply(DLPFCsummary,'model')

makeSortedLogFCplot <- function(x,title,y=NULL){
  #remove duplicate gene identifiers
  x <- x[!duplicated(x$ensembl_gene_id),]
  
  x <- dplyr::mutate(x,direction=(logFC>0))
  x$direction <- gsub(TRUE,'UP',x$direction)
  x$direction <- gsub(FALSE,'DOWN',x$direction)
  if(!is.null(y)){
    x <- dplyr::filter(x,ensembl_gene_id%in%y)
  }
  
  #sort by log fold change
  x <- dplyr::arrange(x,desc(logFC))
  
  
  #fix the factors
  x$ensembl_gene_id <- factor(x$ensembl_gene_id,
                              levels=x$ensembl_gene_id)
  
  g <- ggplot2::ggplot(x,
                       ggplot2::aes(x=ensembl_gene_id,
                                    y=logFC,
                                    color = NULL,
                                    fill = direction,
                                    alpha = significant))
  g <- g+ ggplot2::geom_col()
  g <- g + ggplot2::theme_grey(base_size = 40) 
  g <- g + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
        axis.text.x=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank())
  g <- g + ggplot2::ggtitle(title, subtitle = NULL)

  return(g) 
  #return(x)
}
y <- manifest$aggregateDLPFCblueDLPFC$anno$GeneID
#case and control
png(file='paper_figures/dlpfc_blue_case_control.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(DLPFCsummaryList[[7]],'DLPFC Blue, Case vs Control',y)
dev.off()

png(file='paper_figures/dlpfc_blue_case_control_female.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(DLPFCsummaryList[[11]],'DLPFC Blue, Case vs Control in Female',y)
dev.off()

png(file='paper_figures/dlpfc_blue_case_control_male.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(DLPFCsummaryList[[13]],'DLPFC Blue, Case vs Control in Male',y)
dev.off()

####CER
CBEsummary <- dplyr::select(degResultTable$syn8468023,
                              Model,
                              BrainRegion.ref,
                              Comparison,
                              logFC,
                              adj.P.Val,
                              ensembl_gene_id,
                              hgnc_symbol,
                              Gender.ref,
                              AOD.ref)

TCXsummary <- dplyr::filter(TCXsummary,BrainRegion.ref=='TCX')
TCXsummary$model <- paste0(TCXsummary$Model,'.',TCXsummary$Comparison,'.',TCXsummary$Gender.ref,'.',TCXsummary$AOD.ref)
TCXsummary <- dplyr::select(TCXsummary,-Model,-Comparison,-Gender.ref,-AOD.ref)
TCXsummary <- dplyr::mutate(TCXsummary,significant = adj.P.Val <= 0.05)
TCXsummary <- dplyr::mutate(TCXsummary,gene = ensembl_gene_id)
TCXsummary <- dplyr::select(TCXsummary,-adj.P.Val)

#split into list for each model
TCXsummaryList <- plyr::dlply(TCXsummary,'model')
TCXsummaryList[[8]]$logFC <- -TCXsummaryList[[8]]$logFC

y <- manifest$aggregateTCXblueTCX$anno$GeneID
#case and control
png(file='paper_figures/tcx_blue_case_control.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(TCXsummaryList[[8]],'TCX Blue, Case vs Control',y)
dev.off()

png(file='paper_figures/tcx_blue_case_control_female.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(TCXsummaryList[[10]],'TCX Blue, Case vs Control in Female',y)
dev.off()

png(file='paper_figures/tcx_blue_case_control_apoe.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(TCXsummaryList[[1]],'TCX Blue, APOE E4 Dosage',y)
dev.off()

### cerebellum
CBEsummary <- dplyr::select(degResultTable$syn8468023,
                            Model,
                            BrainRegion.ref,
                            Comparison,
                            logFC,
                            adj.P.Val,
                            ensembl_gene_id,
                            hgnc_symbol,
                            Gender.ref,
                            AOD.ref)

CBEsummary <- dplyr::filter(CBEsummary,BrainRegion.ref=='CER')
CBEsummary$model <- paste0(CBEsummary$Model,'.',CBEsummary$Comparison,'.',CBEsummary$Gender.ref,'.',CBEsummary$AOD.ref)
CBEsummary <- dplyr::select(CBEsummary,-Model,-Comparison,-Gender.ref,-AOD.ref)
CBEsummary <- dplyr::mutate(CBEsummary,significant = adj.P.Val <= 0.05)
CBEsummary <- dplyr::mutate(CBEsummary,gene = ensembl_gene_id)
CBEsummary <- dplyr::select(CBEsummary,-adj.P.Val)

#split into list for each model
CBEsummaryList <- plyr::dlply(CBEsummary,'model')
CBEsummaryList[[8]]$logFC <- -CBEsummaryList[[8]]$logFC

y <- manifest$aggregateCBEturquoiseCBE$anno$GeneID
#case and control
png(file='paper_figures/cbe_blue_case_control.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(CBEsummaryList[[8]],'CBE Turquoise, Case vs Control',y)
dev.off()

png(file='paper_figures/cbe_blue_case_control_female.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(CBEsummaryList[[10]],'CBE Turquoise, Case vs Control in Female',y)
dev.off()

png(file='paper_figures/cbe_blue_case_control_male.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(CBEsummaryList[[11]],'CBE Turquoise, Case vs Control in Male',y)
dev.off()
png(file='paper_figures/cbe_blue_case_control_apoe.png',height=640,width=1280,pointsize=20)
makeSortedLogFCplot(CBEsummaryList[[1]],'CBE Turquoise, APOE E4 Dosage',y)
dev.off()

##### export deg list to file
