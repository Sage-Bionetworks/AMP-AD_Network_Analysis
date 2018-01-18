gef <- read.csv('gef.csv',stringsAsFactors = F)
gtpase <- read.csv('gtpase.csv',stringsAsFactors = F)
synapser::synLogin()
driverScoreObj <- synapser::synGet('syn11688680')
driverScore<-read.csv(driverScoreObj$path,stringsAsFactors = F)
View(driverScore)
View(gef)

gef <- dplyr::left_join(gef,driverScore,by=c('Gene.Name'='external_gene_name'))
write.csv(gef,file='gef_for_paul.csv',quote=F,row.names=F)

gtpase <- dplyr::left_join(gtpase,driverScore,by=c('Gene.Name'='external_gene_name'))
write.csv(gtpase,file='gtpase_for_paul.csv',quote=F,row.names=F)
