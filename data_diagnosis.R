#pull vivek's data
synapseClient::synapseLogin()
foo <- synapseClient::synQuery("select name,id from file where parentId==\'syn10175929\'")
bar <- lapply(foo$file.id,
              synapseClient::synGet)

#pull all mods (ROSMAP)
allMods <- synapseClient::synTableQuery("select * from syn10146524")@values

#pull expression and clinical data (ROSMAP)
library(dplyr)
source('dataPulling/pullExpressionAndPhenoWinsorized.R')

#compute eigengene associations 
listify <- function(x,y,z){
  ###fxn will listify a long form table
  ###x: unique key
  ###y: values
  ###z: keys
  return(unique(y[which(z==x)]))
}
modulesLargeList <- lapply(unique(allMods$ModuleNameFull),
                           listify,
                           allMods$GeneID,
                           allMods$ModuleNameFull)
names(modulesLargeList) <- unique(allMods$ModuleNameFull)



###extract eigengenes
getEigenGenes <- function(geneList,geneExpr){
  w2 <- which(colnames(geneExpr)%in%geneList)
  if(length(w2)>10){
    geneMat <- scale(geneExpr[,w2])
    
    return(svd(geneMat)$u[,1:5])
  }else{
    return(NA)
  }
}

#get ensg genes
geneExpressionForAnalysisEnsg <- lapply(geneExpressionForAnalysis,function(x){
  y <- dplyr::select(x,dplyr::starts_with('ENSG')) %>% data.matrix
  y[is.na(y)] <- 0
  sdvec <- apply(y,2,sd)
  y <- y[,which(sdvec>0)]
  return(y)
})

#test11 <- getEigenGenes(modulesLargeList[[1]],geneExpressionForAnalysisEnsg[[1]])
test11<-utilityFunctions::outerLapply(getEigenGenes,geneExpressionForAnalysisEnsg,modulesLargeList)

adDiagnosis <- lapply(geneExpressionForAnalysis,function(x){
  return(dplyr::select(x,logitDiagnosis))
})

logisticRegressionAggregatePval <- function(x,y){
  #x <- as.numeric(x)
  y <- as.numeric(data.matrix(y))
  #print(x)
  #print(y)
  #print(str(x))
  #print(dim(y))
  logLik1 <- logLik(glm(y~x,family='binomial'))[1]
  #print(logLik1)
  logLik2 <- logLik(glm(y~1,family='binomial'))[1]
  #print(logLik2)
  lrt <- -2*logLik2+2*logLik1
  #print(lrt)
  return(pchisq(lrt,1,lower.tail=F))
}

logisticRegressionOddsRatio <- function(x,y){
  y <- as.numeric(data.matrix(y))
  glmmod <- glm(y~x,family='binomial')
  return(summary(glmmod)$coef[2,1])
}

wrapperFxn1 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn2 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      apply(eigengenes,2,logisticRegressionAggregatePval,diagnosis) %>% return
    }else{
      return(NA)
    }
  }
  lapply(eigengenes,wrapperFxn2,diagnosis) %>% return
}

wrapperFxn3 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn4 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      eigengenes = scale(eigengenes)
      apply(eigengenes,2,logisticRegressionOddsRatio,diagnosis) %>% return
    }else{
      return(NA)
    }
  }
  lapply(eigengenes,wrapperFxn4,diagnosis) %>% return
}

wrapperFxn5 <- function(diagnosis,eigengenes){
  library(dplyr)
  wrapperFxn6 <- function(eigengenes,diagnosis){
    library(dplyr)
    if(!is.na(eigengenes)){
      #eigengenes = scale(eigengenes)
      cor(eigengenes,diagnosis,use = 'pairwise.complete.obs') %>% c %>% return
    }else{
      return(NA)
    }
  }
  res <- lapply(eigengenes,wrapperFxn6,diagnosis)
  names(res) <- names(eigengenes)
  return(res)
}


ModulePvalues <- mapply(wrapperFxn1,
                        adDiagnosis,
                        test11,
                        SIMPLIFY = F)

ModuleOddsRatios <- mapply(wrapperFxn3,
                           adDiagnosis,
                           test11,
                           SIMPLIFY= F)

ModuleCorrelations <- mapply(wrapperFxn5,
                             adDiagnosis,
                             test11,
                             SIMPLIFY = F)
#######focus on rosmap
rosmapEigengenes <- test11$rosmapDLPFC
rosmapExpr <- geneExpressionForAnalysisEnsg$rosmapDLPFC
rosmapPCs <- svd(scale(rosmapExpr))

##### adjust for no pcs
logitWrap <- function(x,y,z=NULL){
  library(dplyr)
  if(!is.null(z)){
    lmres <- glm(y~x+z,family='binomial') %>%
      summary
  }else{
    lmres <- glm(y~x,family='binomial') %>%
      summary    
  }
  return(lmres$coef[2,4])
}
wrapperf <- function(x,y,z=NULL){
  if(!is.na(x)){
    return(apply(x,2,logitWrap,y,z))
  }else{
    return(NA)
  }
}

rosmap_0_pc<-lapply(rosmapEigengenes,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis)
rosmap_0_pc %>%
  unlist %>%
  c %>%
  gap::qqunif()


rosmap_1_pc<-lapply(rosmapEigengenes,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis,rosmapPCs$u[,1])
rosmap_1_pc %>%
  unlist %>%
  c %>%
  gap::qqunif()

rosmap_10_pc<-lapply(rosmapEigengenes,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis,rosmapPCs$u[,1:10])
rosmap_10_pc %>%
  unlist %>%
  c %>%
  gap::qqunif()


rosmapCombined <- dplyr::left_join(geneExpressionForAnalysis$rosmapDLPFC,
                                                          covariateList$syn8456631,
                                                          c('aSampleId'='SampleID'))
rosmapCovariates <- rosmapCombined[,15584:15601]
View(rosmapCovariates)
rosmapCovariates <- dplyr::select(rosmapCovariates,-Diagnosis.x)
rosmapCovariates$Diagnosis.y <- as.factor(rosmapCovariates$Diagnosis.y)
pairs(cbind(rosmapPCs$u[,1:5],rosmapCovariates),horInd=1:5,verInd=6:17)
cor(rosmapPCs$u[,1:5],apply(rosmapCovariates,2,as.numeric),use='pairwise.complete.obs')





covMat <- dplyr::select(rosmapCovariates,-logitDiagnosis,-cogdx,-Diagnosis.y)
covMat$Batch <- as.factor(covMat$Batch)
covMat <- model.matrix(~.-1,covMat)
rosmap_cov<-lapply(rosmapEigengenes,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis,covMat)
rosmap_cov %>%
  unlist %>%
  c %>%
  gap::qqunif(xlim=c(0,3),ylim=c(0,7),col='red')
par(new=T)

rosmap_0_pc %>%
  unlist %>%
  c %>%
  gap::qqunif(xlim=c(0,3),ylim=c(0,7),col='blue')

######## look into eigengenes weights

pcgenemat <- data.frame(GeneID = colnames(rosmapExpr),
                        pc1 = rosmapPCs$v[,1],
                        pc2 = rosmapPCs$v[,2],
                        pc3 = rosmapPCs$v[,3],
                        pc4 = rosmapPCs$v[,4],
                        pc5 = rosmapPCs$v[,5],
                        stringsAsFactors=F)

combpcmods <- dplyr::left_join(allMods,pcgenemat)
genesets1 <- synapseClient::synGet('syn5923958')
load(synapseClient::getFileLocation(genesets1))
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]

cellList <- GeneSets$Cell_Markers
cellList$ADGWAS <- adList
combpcmods$mg <- as.factor(combpcmods$external_gene_name%in%cellList$`Zhang:Microglia`)
combpcmods$endo <- as.factor(combpcmods$external_gene_name%in%cellList$`Zhang:Endothelial`)
combpcmods$astro <- as.factor(combpcmods$external_gene_name%in%cellList$`Zhang:Astrocyte`)
combpcmods$neuron <- as.factor(combpcmods$external_gene_name%in%cellList$`Zhang:Neuron`)
combpcmods$myeloligo <- as.factor(combpcmods$external_gene_name%in%cellList$`Zhang:MyelinOligos`)
pairs(combpcmods[,c('mg','endo','astro','neuron','myeloligo','pc1','pc2','pc3','pc4','pc5')],horInd = 1:5,verInd = 6:10)
cellType <- as.numeric(combpcmods$mg)
cellType[which(cellType==2)] <- 'blue'
cellType[which(as.numeric(combpcmods$endo)==2)] <- 'red'
cellType[which(as.numeric(combpcmods$astro)==2)] <- 'green'
cellType[which(as.numeric(combpcmods$neuron)==2)] <- 'purple'
cellType[which(as.numeric(combpcmods$myeloligo)==2)] <- 'orange'
cellType[which(cellType==1)] <- NA
pairs(combpcmods[,c('pc2','pc3','pc4')],
      col=(cellType),
      pch=15)

####pull rwgcna modules
rwgcna <- synapseClient::synTableQuery("SELECT * FROM syn10163855 WHERE ( ( brainRegion = 'DLPFC' ) AND ( method = 'rWGCNA' ) )")@values

View(rwgcna)
rwgcnaLargeList <- lapply(unique(rwgcna$ModuleNameFull),
                           listify,
                           rwgcna$GeneID,
                           rwgcna$ModuleNameFull)
names(rwgcnaLargeList) <- unique(rwgcna$ModuleNameFull)
test22<-utilityFunctions::outerLapply(getEigenGenes,geneExpressionForAnalysisEnsg,rwgcnaLargeList)
rosmapEigengenesRwgcna <- test22$rosmapDLPFC

rosmap_0_pc_rwgcna<-lapply(rosmapEigengenesRwgcna,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis)
rosmap_0_pc_rwgcna %>%
  unlist %>%
  c %>%
  gap::qqunif()


rosmap_5_pc_rwgcna<-lapply(rosmapEigengenesRwgcna,wrapperf,adDiagnosis$rosmapDLPFC$logitDiagnosis,rosmapPCs$u[,1:5])
rosmap_5_pc_rwgcna %>%
  unlist %>%
  c %>%
  gap::qqunif()

#compare to associations in data that was provided by Vivek



