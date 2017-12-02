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

set.seed(1)
agg_pval[agg_pval==1] <- runif(sum(agg_pval==1))
ind_pval[ind_pval==1] <- runif(sum(ind_pval==1))
set.seed(1)
foobar2$pval[foobar2$pval==1] <- runif(sum(foobar2$pval==1))

adGeneticsSummaryAgg$GeneSetAssociationStatistic <- agg_pval
adGeneticsSummaryInd$GeneSetAssociationStatistic <- ind_pval


gap::qqunif(dplyr::filter(adGeneticsSummaryInd,GeneSetName=='genecards')$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42))
par(new=T)
gap::qqunif(dplyr::filter(foobar2,pathway=='genecards')$pval,xlim=c(0,5),ylim=c(0,42),col='green')
par(new=T)
gap::qqunif(dplyr::filter(adGeneticsSummaryAgg,GeneSetName=='genecards')$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red')

gap::qqunif(adGeneticsSummaryInd$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42))
par(new=T)
gap::qqunif(foobar2$pval,xlim=c(0,5),ylim=c(0,42),col='green')
par(new=T)
gap::qqunif(adGeneticsSummaryAgg$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red',main='AD Gene-set Enrichments')
legend('topleft',c('all modules','all DEGs','aggregate modules'),fill=c('blue','green','red'))

ks.test(adGeneticsSummaryAgg$GeneSetAssociationStatistic,
        foobar2$pval,alternative='greater')
