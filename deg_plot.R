#make volcano plot for talk
synapseClient::synapseLogin()

foo <- synapseClient::synGet('syn8456721')
bar <- data.table::fread(foo@filePath,data.table=F)
bar1 <- dplyr::filter(bar,Model=='Diagnosis')
bar1 <- dplyr::filter(bar1,Comparison=='AD-CONTROL')

p = ggplot2::ggplot(bar1, 
                    ggplot2::aes(y = -log10(adj.P.Val), 
                                 x = logFC, color = Direction)) + ggplot2::geom_point() + ggplot2::xlim(c(-1,1))
p = p + ggplot2::scale_color_manual(values = c('green','grey','red'))
#p = p + ggplot2::facet_grid(.~Comparison, scales = 'fixed')
p
