synapseClient::synapseLogin()
foo <- synapseClient::synTableQuery("SELECT * FROM syn9770842")@values


foo1<-tidyr::gather(foo,test,pval,eigengene1,eigengene2,eigengene3,eigengene4,eigengene5,eigengeneAggregate,adGeneticEnrich,microgEnrich,neurEnrich,astroEnrich,endoEnrich)
a4<-(which(!duplicated(dplyr::select(foo1,-brainRegionAssociation))))
foo1 <- foo1[a4,]
foo1b <- dplyr::filter(foo1,test=='adGeneticEnrich')

foo1b <- dplyr::filter(foo1b,p.adjust(pval,method='fdr')<=0.05)

bb<-table(foo1b$method)
aa<-c('speakEasy'=280,
      'megena'=3155,
      'metanetwork'=96,
      'wina'=594,
      'consensus'=463)
cc <- cbind(aa[names(bb)],bb)
par(mar=c(6,4,4,2))
barplot((cc[,2]/cc[,1])*100,las=2,ylab='percent modules enriched',main='% modules enriched for AD Genetic Loci')
table(foo1b$brainRegion)

consensusModule60 <- dplyr::filter(foo1,ModuleNameFull=='consensus60CER')
par(mar=c(8,6,4,4))
barplot(-log10(consensusModule60$pval)[36:47],
        las=2,
        log='y',
        names.arg = c(consensusModule60$brainRegionAssociation[36:42],consensusModule60$test[43:47]),
        ylab='-log10 p-value',
        main = 'module 60, consensus, CER')
lines(x=c(0,20),y=-log10(c(5e-4,5e-4)),col='red',lwd=2)
