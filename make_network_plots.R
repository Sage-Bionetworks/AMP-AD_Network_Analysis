load('aggregate_module_mainfest.rda')


manifest <- fullManifest

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
filenames <- paste0('./edgeListFolder/',names(edgeLists),'.csv')
mapply(function(x,y){
  write.csv(x[,3:4],file=y,quote=F)
},edgeLists,filenames)

synapseClient::synapseLogin()
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
cellTypeDf <- utilityFunctions::list2df(GeneSets$Cell_Markers)
write.csv(cellTypeDf,file='cellTypeAnnotations.csv',quote=F)


pullReferenceGeneSets <- function(url){
  #pull gene set from Enrichr
  foo <- RCurl::getURL(url)
  
  #split by carriage returns
  fooSplit <- strsplit(foo,"\n")
  
  #unlist
  fooSplit <- unlist(fooSplit)
  
  #split by tab
  fooSplit <- strsplit(fooSplit,"\t")
  
  #name them by the first element
  catNames <- lapply(fooSplit,function(x) x[1])
  names(fooSplit) <- catNames
  
  #drop the first two elements
  fooSplit <- lapply(fooSplit,function(x) x[-c(1,2)])
  
  fxn1 <- function(x){
    library(dplyr)
    strsplit(x,',') %>%
      sapply(function(y) y[1]) %>%
      return
  }
  #remove 1.0 from elements
  fooSplit2 <- lapply(fooSplit,fxn1)
  names(fooSplit2) <- names(fooSplit)
  return(fooSplit2)
}
#splicing and mitochondrial dna repair
#mRNA.splicing..via.spliceosome
go_bp_url <- "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2017b"
go_bp <- pullReferenceGeneSets(go_bp_url)
panther <- pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")

go_keep <- go_bp[c('mRNA splicing, via spliceosome','mitochondrial DNA repair','detection of unfolded protein','histone H3-K9 dimethylation')]
go_keep$mitochondrion <- go_cc$mitochondrion
go_keep$ubiquitin <- panther$`Ubiquitin proteasome pathway_Homo sapiens_P00060`
#go_bp$`mRNA splicing, via spliceosome`
names(go_keep)[1:4] <- c('splicing','mtdnarepair','unfolded','histonh3k9')
annoDf <- utilityFunctions::list2df(go_keep)
write.csv(annoDf,file='splicing_mtdna.csv',quote=F)


######load degs

library(dplyr)
synapseClient::synapseLogin()
degResProv<-synapseClient::synGetActivity('syn10496554')
synapseIdsOfDegs <- sapply(degResProv@properties$used,function(x) x$reference$targetId) %>% unlist

synapseIdsOfDegs[synapseIdsOfDegs == 'syn10157628'] <- 'syn10526259'
degResultTable <- rSynapseUtilities::loadDelimIntoList(synapseIdsOfDegs)


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

TCXsummary <- dplyr::select(degResultTable$syn8468023,
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

PHGsummary <- dplyr::select(degResultTable$syn10526259,
                            Model,
                            BrainRegion,
                            Comparison,
                            logFC,
                            adj.P.Val,
                            ensembl_gene_id,
                            hgnc_symbol,
                            reference,
                            against)

#PHGsummary <- dplyr::filter(PHGsummary,BrainRegion=='PHG')
wphg<-union(grep('PHG',PHGsummary$Comparison),grep('PHG',PHGsummary$BrainRegion))
PHGsummary <- PHGsummary[wphg,]
PHGsummary$model <- paste0(PHGsummary$Model,'.',PHGsummary$Comparison,'.',PHGsummary$reference,'.',PHGsummary$against)
PHGsummary <- dplyr::select(PHGsummary,-Model,-Comparison,-reference,-against)
PHGsummary <- dplyr::mutate(PHGsummary,significant = adj.P.Val <= 0.05)
PHGsummary <- dplyr::mutate(PHGsummary,gene = ensembl_gene_id)
PHGsummary <- dplyr::select(PHGsummary,-adj.P.Val)

PHGsummaryList <- plyr::dlply(PHGsummary,'model')
PHGsummaryList[[8]]$logFC <- -PHGsummaryList[[8]]$logFC


####STG
STGsummary <- dplyr::select(degResultTable$syn10526259,
                            Model,
                            BrainRegion,
                            Comparison,
                            logFC,
                            adj.P.Val,
                            ensembl_gene_id,
                            hgnc_symbol,
                            reference,
                            against)

#STGsummary <- dplyr::filter(STGsummary,BrainRegion=='STG')
wSTG<-union(grep('STG',STGsummary$Comparison),grep('STG',STGsummary$BrainRegion))
STGsummary <- STGsummary[wSTG,]
STGsummary$model <- paste0(STGsummary$Model,'.',STGsummary$Comparison,'.',STGsummary$reference,'.',STGsummary$against)
STGsummary <- dplyr::select(STGsummary,-Model,-Comparison,-reference,-against)
STGsummary <- dplyr::mutate(STGsummary,significant = adj.P.Val <= 0.05)
STGsummary <- dplyr::mutate(STGsummary,gene = ensembl_gene_id)
STGsummary <- dplyr::select(STGsummary,-adj.P.Val)

STGsummaryList <- plyr::dlply(STGsummary,'model')
STGsummaryList[[8]]$logFC <- -STGsummaryList[[8]]$logFC


###IFG
IFGsummary <- dplyr::select(degResultTable$syn10526259,
                            Model,
                            BrainRegion,
                            Comparison,
                            logFC,
                            adj.P.Val,
                            ensembl_gene_id,
                            hgnc_symbol,
                            reference,
                            against)

#IFGsummary <- dplyr::filter(IFGsummary,BrainRegion=='IFG')
wIFG<-union(grep('IFG',IFGsummary$Comparison),grep('IFG',IFGsummary$BrainRegion))
IFGsummary <- IFGsummary[wIFG,]
IFGsummary$model <- paste0(IFGsummary$Model,'.',IFGsummary$Comparison,'.',IFGsummary$reference,'.',IFGsummary$against)
IFGsummary <- dplyr::select(IFGsummary,-Model,-Comparison,-reference,-against)
IFGsummary <- dplyr::mutate(IFGsummary,significant = adj.P.Val <= 0.05)
IFGsummary <- dplyr::mutate(IFGsummary,gene = ensembl_gene_id)
IFGsummary <- dplyr::select(IFGsummary,-adj.P.Val)

IFGsummaryList <- plyr::dlply(IFGsummary,'model')
IFGsummaryList[[8]]$logFC <- -IFGsummaryList[[8]]$logFC


####FP
FPsummary <- dplyr::select(degResultTable$syn10526259,
                            Model,
                            BrainRegion,
                            Comparison,
                            logFC,
                            adj.P.Val,
                            ensembl_gene_id,
                            hgnc_symbol,
                            reference,
                            against)

#FPsummary <- dplyr::filter(FPsummary,BrainRegion=='FP')
wFP<-union(grep('FP',FPsummary$Comparison),grep('FP',FPsummary$BrainRegion))
FPsummary <- FPsummary[wFP,]
FPsummary$model <- paste0(FPsummary$Model,'.',FPsummary$Comparison,'.',FPsummary$reference,'.',FPsummary$against)
FPsummary <- dplyr::select(FPsummary,-Model,-Comparison,-reference,-against)
FPsummary <- dplyr::mutate(FPsummary,significant = adj.P.Val <= 0.05)
FPsummary <- dplyr::mutate(FPsummary,gene = ensembl_gene_id)
FPsummary <- dplyr::select(FPsummary,-adj.P.Val)

FPsummaryList <- plyr::dlply(FPsummary,'model')
FPsummaryList[[8]]$logFC <- -FPsummaryList[[8]]$logFC


###write out case v control logfc
write.csv(DLPFCsummaryList[[10]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='dlpfc_case_v_control.csv',quote=F)

write.csv(TCXsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='tcx_case_v_control.csv',quote=F)

write.csv(CBEsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='cbe_case_v_control.csv',quote=F)

write.csv(FPsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='fp_case_v_control.csv',quote=F)

write.csv(STGsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='stg_case_v_control.csv',quote=F)

write.csv(IFGsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='ifg_case_v_control.csv',quote=F)

write.csv(PHGsummaryList[[8]][,c('logFC','ensembl_gene_id','hgnc_symbol','significant')],file='phg_case_v_control.csv',quote=F)


####write out female case vs control log fc
foo <- DLPFCsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='dlpfc_case_v_control_female.csv',quote=F)

foo <- TCXsummaryList[[10]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='tcx_case_v_control_female.csv',quote=F)

foo <- CBEsummaryList[[10]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='cbe_case_v_control_female.csv',quote=F)

foo <- FPsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='fp_case_v_control_female.csv',quote=F)

foo <- STGsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='stg_case_v_control_female.csv',quote=F)

foo <- IFGsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='ifg_case_v_control_female.csv',quote=F)

foo <- PHGsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_female','hgnc_symbol','significant_female')
write.csv(foo,file='phg_case_v_control_female.csv',quote=F)


####write out male case vs control log fc
foo <- DLPFCsummaryList[[13]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='dlpfc_case_v_control_male.csv',quote=F)

foo <- TCXsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='tcx_case_v_control_male.csv',quote=F)

foo <- CBEsummaryList[[11]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='cbe_case_v_control_male.csv',quote=F)

foo <- FPsummaryList[[12]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='fp_case_v_control_male.csv',quote=F)

foo <- STGsummaryList[[12]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='stg_case_v_control_male.csv',quote=F)

foo <- IFGsummaryList[[12]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='ifg_case_v_control_male.csv',quote=F)

foo <- PHGsummaryList[[12]][,c('logFC','hgnc_symbol','significant')]
colnames(foo) <- c('logFC_male','hgnc_symbol','significant_male')
write.csv(foo,file='phg_case_v_control_male.csv',quote=F)

