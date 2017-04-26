synapseClient::synapseLogin()

#get synIds for gene expression variables
geneExpressionDataManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( dataSubType = 'residualGeneExpForNetAnlz' ) AND ( normalizationType = 'CQN' ) )")

#get synIds for covariates
covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

#geneExpressionDataObj <- sapply(geneExpressionDataManifest@values$id,synapseClient::synGet)
#covariateManifestObj <- sapply(covariateManifest@values$id,synapseClient::synGet)

#load expression data into R
geneExpressionList <- rSynapseUtilities::loadDelimIntoList(geneExpressionDataManifest@values$id)

#load covariate data into R
covariateList <- rSynapseUtilities::loadDelimIntoList(covariateManifest@values$id)


#split mayo into two data-frames
geneExpressionForAnalysis <- list()
geneExpressionForAnalysis$mayoTCX <- dplyr::select(geneExpressionList$syn8466826,
                         dplyr::ends_with('TCX'))
rownames(geneExpressionForAnalysis$mayoTCX) <- geneExpressionList$syn8466826$ensembl_gene_id

geneExpressionForAnalysis$mayoCER <- dplyr::select(geneExpressionList$syn8466826,
                         dplyr::ends_with('CER'))
rownames(geneExpressionForAnalysis$mayoCER) <- geneExpressionList$syn8466826$ensembl_gene_id
library(dplyr)
#rosmap
geneExpressionForAnalysis$rosmapDLPFC <- geneExpressionList$syn8456719 %>% dplyr::select(-ensembl_gene_id)
rownames(geneExpressionForAnalysis$rosmapDLPFC) <- geneExpressionList$syn8456719$ensembl_gene_id
#split mssm into 4 data-frames
mssmcovObj <- synapseClient::synGet('syn6100548')
mssmcov <- data.table::fread(mssmcovObj@filePath,data.table=F)
fpIds <- dplyr::filter(mssmcov,BrodmannArea=='BM10')%>%
  dplyr::select(sampleIdentifier)
stgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM22')%>%
  dplyr::select(sampleIdentifier)
phgIds <- dplyr::filter(mssmcov,BrodmannArea=='BM36')%>%
  dplyr::select(sampleIdentifier)
IFGIds <- dplyr::filter(mssmcov,BrodmannArea=='BM44')%>%
  dplyr::select(sampleIdentifier)

geneExpressionForAnalysis$msbbFP <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%fpIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbSTG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%stgIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbPHG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%phgIds$sampleIdentifier)]
geneExpressionForAnalysis$msbbIFG <- geneExpressionList$syn8485027[,which(colnames(geneExpressionList$syn8485027)%in%IFGIds$sampleIdentifier)]
rownames(geneExpressionForAnalysis$msbbFP) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbSTG) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbPHG) <- geneExpressionList$syn8485027$ensembl_gene_id
rownames(geneExpressionForAnalysis$msbbIFG) <- geneExpressionList$syn8485027$ensembl_gene_id

#explode by tissue