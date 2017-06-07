synapseClient::synapseLogin()

#pull DEGs
foo <- synapseClient::synGet('syn9766351')
table(all.results$CER$Comparison)

#restrict to just AD-CONTROLs
classic <- lapply(all.results,function(x){
  return(dplyr::filter(x,Comparison=='AD-CONTROL'))
})

#get positive with abs(log(FC))>1.2
positive <- lapply(classic,function(x){
  return(dplyr::filter(x,p.adjust(P.Value,method='fdr')<=0.05))
})

negative <- lapply(classic,function(x){
  return(dplyr::filter(x,p.adjust(P.Value,method='fdr')<=0.05))
})

