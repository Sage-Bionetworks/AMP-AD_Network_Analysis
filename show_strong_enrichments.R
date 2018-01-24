synapseClient::synapseLogin()

source('moduleAnalysis/summaryManifestFunctions.R')

adGeneticsSummaryAgg <- getAdGenetics(synId='syn10915669')
adGeneticsSummaryInd <- getAdGenetics(synId='syn10309369')
adGeneticsSummaryAgg2 <- getAdGenetics2(synId='syn10915669')

adList<-adGeneticsSummaryAgg2[[2]]

degResObj <- synapseClient::synGet("syn10496554")
load(degResObj@filePath)

foo2 <- synapseClient::synTableQuery("select distinct external_gene_name from syn10309369")@values

foobar <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperPval, amp.ad.de.geneSets, adList,foo2$external_gene_name)

foobar <- data.frame(foobar,stringsAsFactors=F)
foobar$pathway <- rownames(foobar)
foobar2<-tidyr::gather(foobar,key='geneset',value='pval',-pathway)

foobar3 <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperOR, amp.ad.de.geneSets, adList,foo2$external_gene_name)
foobar3 <- data.frame(foobar3,stringsAsFactors=F)
foobar3$pathway <- rownames(foobar)
foobar4<-tidyr::gather(foobar3,key='geneset',value='OR',-pathway)

agg_eff <- adGeneticsSummaryAgg$GeneSetEffect
ind_eff <- adGeneticsSummaryInd$GeneSetEffect
de_eff <- foobar4$OR

agg_eff[!is.finite(agg_eff)] <- NA
ind_eff[!is.finite(ind_eff)] <- NA
de_eff[!is.finite(de_eff)] <- NA
agg_eff[agg_eff==0] <- NA
ind_eff[ind_eff==0] <- NA
de_eff[de_eff==0] <- NA

hist(log10(adGeneticsSummaryInd$GeneSetEffect))
mean(log10(adGeneticsSummaryInd$GeneSetEffect))


agg_pval <- adGeneticsSummaryAgg$GeneSetAssociationStatistic
ind_pval <- adGeneticsSummaryInd$GeneSetAssociationStatistic

# set.seed(1)
# agg_pval[agg_pval==1] <- runif(sum(agg_pval==1))
# ind_pval[ind_pval==1] <- runif(sum(ind_pval==1))
# set.seed(1)
# foobar2$pval[foobar2$pval==1] <- runif(sum(foobar2$pval==1))
# 
# adGeneticsSummaryAgg$GeneSetAssociationStatistic <- agg_pval
# adGeneticsSummaryInd$GeneSetAssociationStatistic <- ind_pval

gaiteri_mods <- read.csv('zhang_modules.csv',stringsAsFactors=F)
modList <- lapply(unique(gaiteri_mods$Module), utilityFunctions::listify,gaiteri_mods$Gene_Symbol,gaiteri_mods$Module)
names(modList) <- unique(gaiteri_mods$Module)
foobarX <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperPval, modList, adList,foo2$external_gene_name)

foobarX <- data.frame(foobarX,stringsAsFactors=F)
foobarX$pathway <- rownames(foobarX)
foobarX2<-tidyr::gather(foobarX,key='geneset',value='pval',-pathway)

foobarX3 <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperOR, modList, adList,foo2$external_gene_name)
foobarX3 <- data.frame(foobarX3,stringsAsFactors=F)
foobarX3$pathway <- rownames(foobarX)
foobarX4<-tidyr::gather(foobarX3,key='geneset',value='OR',-pathway)






gap::qqunif(dplyr::filter(adGeneticsSummaryInd,GeneSetName=='genecards')$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42))
par(new=T)
gap::qqunif(dplyr::filter(foobar2,pathway=='genecards')$pval,xlim=c(0,5),ylim=c(0,42),col='green')
par(new=T)
gap::qqunif(dplyr::filter(adGeneticsSummaryAgg,GeneSetName=='genecards')$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red')



uniqueMods <- c(CBEres$moduleGraph$from,CBEres$moduleGraph$to)%>% unique
uniqueMods <- c(DLPFCres$moduleGraph$from,DLPFCres$moduleGraph$to) %>% unique %>% c(uniqueMods)
uniqueMods <- c(PHGres$moduleGraph$from,PHGres$moduleGraph$to) %>% unique %>% c(uniqueMods)
uniqueMods <- c(FPres$moduleGraph$from,FPres$moduleGraph$to) %>% unique %>% c(uniqueMods)
uniqueMods <- c(TCXres$moduleGraph$from,TCXres$moduleGraph$to) %>% unique %>% c(uniqueMods)
uniqueMods <- c(IFGres$moduleGraph$from,IFGres$moduleGraph$to) %>% unique %>% c(uniqueMods)
uniqueMods <- c(STGres$moduleGraph$from,STGres$moduleGraph$to) %>% unique %>% c(uniqueMods)


gap::qqunif(adGeneticsSummaryInd$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42))
par(new=T)
gap::qqunif(dplyr::filter(adGeneticsSummaryInd,ModuleNameFull %in% uniqueMods)$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='purple')
par(new=T)
gap::qqunif(foobar2$pval,xlim=c(0,5),ylim=c(0,42),col='green')
par(new=T)
gap::qqunif(foobarX2$pval,xlim=c(0,5),ylim=c(0,42),col='brown')
par(new=T)
gap::qqunif(adGeneticsSummaryAgg$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red',main='AD Gene-set Enrichments')
legend('topleft',c('all modules','all DEG enriched modules','all DEGs','aggregate modules'),fill=c('blue','purple','green','red'))

ks.test(adGeneticsSummaryAgg$GeneSetAssociationStatistic,
        foobar2$pval,alternative='greater')
