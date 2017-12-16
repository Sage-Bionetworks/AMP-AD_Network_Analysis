#pull trena results from synapse

synapser::synLogin()
trena_results <- synapser::synGet('syn11376251')
foo <- load(trena_results$path)

#filter mods to be brain region specific
mods <- list()
mods$cer_mods <- sage10.mayo.cer[grep('CBE',names(sage10.mayo.cer))]
mods$tcx_mods <- sage10.mayo.tcx[grep('TCX',names(sage10.mayo.tcx))]
mods$fp_mods <- sage10.mssm.fp[grep('FP',names(sage10.mssm.fp))]
mods$ifg_mods <- sage10.mssm.ifg[grep('IFG',names(sage10.mssm.ifg))]
mods$phg_mods <- sage10.mssm.phg[grep('PHG',names(sage10.mssm.phg))]
mods$stg_mods <- sage10.mssm.phg[grep('STG',names(sage10.mssm.stg))]
mods$dlpfc_mods <- sage10.rosmap[grep('DLPFC',names(sage10.rosmap))]

mods <- Reduce(c,mods)


#extract top 10 driver genes for each module
get_top_10 <- function(x){
  return(x$gene[1:10])
}

geneLists<-lapply(mods,get_top_10)
geneListDf <- utilityFunctions::list2df(geneLists)
write.csv(geneListDf,file='trenaTopDrivers.csv',quote=F,row.names=F)
View(geneListDf)

permLink =githubr::getPermlink(repository = 'Sage-Bionetworks/AMP-AD_Network_Analysis',
                               ref = 'branch',
                               refName = 'module_ranking',
                               repositoryPath = 'get_trena_results.R')

#store to synapse
bar<-synapser::File('trenaTopDrivers.csv',
               parentId='syn11376150',
               versionComment='top 10 trena drivers for each aggregate module')
bar <- synapser::synStore(bar, 
                   used = 'syn11376251',
                   executed = permLink)


