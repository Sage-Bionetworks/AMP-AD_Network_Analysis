#login to Synapse
synapseClient::synapseLogin()

#pull enrichments
enrichments <- synapseClient::synTableQuery("SELECT * FROM syn10492048")@values
enrichments2 <- dplyr::select(enrichments,ModuleNameFull,category,geneSet,fisherPval,fisherOR)
colnames(enrichments2) <- c('ModuleNameFull',
                            'GeneSetName',
                            'GeneSetCategoryName',
                            'GeneSetAssociationStatistic',
                            'GeneSetEffect')

enrichments2$GeneSetBrainRegion <- rep(NA,nrow(enrichments2))
enrichments2$GeneSetDirectionAD <- rep(NA,nrow(enrichments2))
ad_lists <- grep('alzheimer',(enrichments2$GeneSetName),ignore.case=TRUE)
#ad_lists2 <- grep('load',unique(enrichments2$GeneSetName))
#ad_lists3 <- grep("AD",unique(enrichments2$GeneSetName))
#ad_lists<-grep('',enrichments2$GeneSetName)
enrichments2$GeneSetADLinked <- rep(FALSE,nrow(enrichments2))
enrichments2$GeneSetADLinked[ad_lists] <- TRUE

#bonferroni womp womp

bonferroni_fun <- function(x,ntests=1e8){
  return(min(1,x*ntests))
}

enrichments2$GeneSetAdjustedAssociationStatistic <- sapply(enrichments2$GeneSetAssociationStatistic,bonferroni_fun)
enrichments3 <- dplyr::filter(enrichments2,GeneSetAdjustedAssociationStatistic <= 0.05)
View(enrichments3)

boo3 <- dplyr::group_by(enrichments3,ModuleNameFull,GeneSetCategoryName)
boo4 <- dplyr::summarise(boo3,mean_or=mean(GeneSetEffect),mean_pval=mean(-log10(GeneSetAssociationStatistic)),nsigcat=length(GeneSetCategoryName))

### get module sizes
allMods <- synapseClient::synTableQuery("SELECT * FROM syn10338156")@values
library(dplyr)
fooSummarize <- dplyr::group_by(allMods,brainRegion,ModuleName,method)%>%
  dplyr::summarise(numberOfGenes=length(ModuleName))

fooSummarize <- dplyr::mutate(fooSummarize,ModuleNameFull = paste0(ModuleName,brainRegion))
module_score_manifest <- dplyr::left_join(boo4,fooSummarize)

module_score_manifest$mean_or[!is.finite(module_score_manifest$mean_or)] <- NA
module_score_manifest<-na.omit(module_score_manifest)
module_score_manifest$adjusted_log10_mean_or <- MASS::rlm(log10(module_score_manifest$mean_or) ~ log10(module_score_manifest$numberOfGenes) -1)$resid
module_score_manifest$adjusted_log10_nsigset <- MASS::rlm(log10(module_score_manifest$nsigcat) ~ log10(module_score_manifest$numberOfGenes) - 1)$resid
module_score_manifest$adjusted_mean_pval <- MASS::rlm(module_score_manifest$mean_pval ~ log10(module_score_manifest$numberOfGenes) - 1)$resid

module_score_manifest$adjusted_mean_or <- 10^(MASS::rlm(log10(module_score_manifest$mean_or) ~ log10(module_score_manifest$numberOfGenes) -1)$resid)
module_score_manifest$adjusted_nsigset <- 10^(MASS::rlm(log10(module_score_manifest$nsigcat) ~ log10(module_score_manifest$numberOfGenes) - 1)$resid)

module_score_manifest <- dplyr::mutate(module_score_manifest,log10_mean_or = log10(mean_or))
module_score_manifest <- dplyr::mutate(module_score_manifest,log10_nsigcat = log10(nsigcat))
#####box plots with dots
#par(mfcol=c(2,1))
makeBoxPlot <- function(msm,
                        yaxis,
                        yaxisLabel){
  #str <- paste0("y=log10(",yaxis,")")
  g <- ggplot2::ggplot(msm, 
                       ggplot2::aes_string(x="method",
                                    y=yaxis,
                                    color="method",
                                    fill="method"))
  g <- g + ggplot2::geom_boxplot(alpha = 0.2)
  g <- g + ggplot2::geom_jitter(ggplot2::aes_string(x="method",
                                             y=yaxis,
                                             color="method",
                                             fill="method"))
  #g <- g + ggplot2::scale_y_log10()
  g <- g + ggplot2::theme_grey(base_size = 20) 
  g <- g + ggplot2::coord_flip()

  return(g)
}
makeBoxPlot(module_score_manifest,"log10_mean_or")
makeBoxPlot(module_score_manifest,"adjusted_log10_mean_or")

makeBoxPlot(module_score_manifest,"log10_nsigcat")
makeBoxPlot(module_score_manifest,"adjusted_log10_nsigset")

makeBoxPlot(module_score_manifest,"mean_pval")
makeBoxPlot(module_score_manifest,"adjusted_mean_pval")

#####fix msm
msm_new <- dplyr::select(module_score_manifest,
                         ModuleNameFull,
                         method,
                         GeneSetCategoryName,
                         log10_mean_or,
                         mean_pval,
                         log10_nsigcat,
                         adjusted_log10_mean_or,
                         adjusted_mean_pval,
                         adjusted_log10_nsigset)

msm_new2 <- tidyr::gather(msm_new,
                          key_or,
                          value_or,
                          log10_mean_or,
                          adjusted_log10_mean_or)

msm_new2 <- tidyr::gather(msm_new2,
                          key_pval,
                          value_pval,
                          mean_pval,
                          adjusted_mean_pval)

msm_new2 <- tidyr::gather(msm_new2,
                          key_cat,
                          value_cat,
                          log10_nsigcat,
                          adjusted_log10_nsigset)
makeBoxPlot <- function(msm,
                        yaxis,
                        key,
                        category){
  g <- ggplot2::ggplot(dplyr::filter(msm_new2,GeneSetCategoryName==category), 
                       ggplot2::aes_string(x="method",
                                           y=yaxis,
                                           color=key,
                                           fill=key))
  g <- g + ggplot2::geom_boxplot(alpha = 0.2,
                                 position = 'dodge')
  # g <- g + ggplot2::geom_jitter(ggplot2::aes_string(x="method",
  #                                                    y="value_cat",
  #                                                    color="key_cat",
  #                                                    fill="key_cat"))
  
  #g <- g + ggplot2::scale_y_log10()
  g <- g + ggplot2::theme_grey(base_size = 20) 
  #g <- g + ggplot2::coord_flip()
  g
}
makeBoxPlot(msm_new2,"value_cat","key_cat","panther")
