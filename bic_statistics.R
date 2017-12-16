synapseClient::synapseLogin()
foobar <- synapseClient::synGet('syn8268669')
load(foobar@filePath)
fob <- as.matrix(bicNetworks$network)
fob2 <- igraph::graph_from_adjacency_matrix(fob,mode='undirected')
fob3 <- igraph::degree(fob2)

library(synapser)
synapser::synLogin()
genesets1 <- synapser::synGet('syn5923958')
load(genesets1$path)
adTypeDf <- utilityFunctions::list2df(GeneSets$Alzheimers)

mapFun <- utilityFunctions::convertEnsemblToHgnc(names(fob3))
mapFun2 <- utilityFunctions::convertHgncToEnsembl(adTypeDf$value)

colnames(mapFun)
colnames(adTypeDf)

adTypeDf <- dplyr::left_join(adTypeDf,mapFun2,by=c('value'='external_gene_name'))
View(adTypeDf)
adTypeDf <- adTypeDf[!duplicated(adTypeDf),]
View(adTypeDf)

hubDf <- data.frame(ensembl_gene_id=names(fob3),
                    degree=fob3,
                    stringsAsFactors = F)

adTypeDf <- dplyr::left_join(adTypeDf,hubDf)
View(adTypeDf)

keep <- c('MouseMicroglia:2month_TG_vs_WT',
          'MouseMicroglia:4month_TG_vs_WT',
          'MouseMicroglia:6month_TG_vs_WT',
          'MouseMicroglia:8month_TG_vs_WT')

redAdTypeDf <- dplyr::filter(adTypeDf,key%in%keep)

gpl <- ggplot2::ggplot(redAdTypeDf,
                       ggplot2::aes(x=key,
                                    y=degree))
gpl <- gpl + ggplot2::geom_boxplot()
gpl <- gpl + ggplot2::geom_jitter(alpha=0.2)
gpl <- gpl + ggplot2::scale_x_discrete(name ="rTg4501 DEG set", 
                              labels=c('MouseMicroglia:2month_TG_vs_WT'="2 Months",
                                       'MouseMicroglia:4month_TG_vs_WT'="4 Months",
                                       'MouseMicroglia:6month_TG_vs_WT'="6 Months",
                                       'MouseMicroglia:8month_TG_vs_WT'="8 Months"))
gpl <- gpl + ggplot2::scale_y_continuous(name = 'Degree of Gene in Network')
gpl
summary(lm(degree ~ key,data=redAdTypeDf))
mo2 <- dplyr::filter(redAdTypeDf,key=='MouseMicroglia:2month_TG_vs_WT')$degree
mo4 <-dplyr::filter(redAdTypeDf,key=='MouseMicroglia:4month_TG_vs_WT')$degree
mo6 <-dplyr::filter(redAdTypeDf,key=='MouseMicroglia:6month_TG_vs_WT')$degree
mo8 <-dplyr::filter(redAdTypeDf,key=='MouseMicroglia:8month_TG_vs_WT')$degree
wilcox.test(mo2,c(mo4,mo6,mo8),alternative = 'greater')
