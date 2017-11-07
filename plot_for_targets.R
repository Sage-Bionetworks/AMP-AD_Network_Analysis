foobar <- data.table::fread('targetDf.csv',data.table=F)


c2n <- function(y){
  convert2numeric <- function(x){
    if(x=='low'){
      return(0)
    }else if(x=='medium'){
      return(1)
    }else if(x=='high'){
      return(2)
    }
  }
  return(sapply(y,convert2numeric))
}
rownames(foobar) <- foobar$NominatedTarget
foobar <- foobar[,-1]
foobar[1:5,]
foobar <- apply(foobar,2,c2n)

pheatmap::pheatmap(foobar)

svd2 <- svd(scale(foobar))
#pairs(svd2$u[,1:3])
cross_consortia_score <- rowSums(foobar)
sort(cross_consortia_score,decreasing=T)

foobar2 <- read.csv('moleculaEvidenceClassSummary.csv',stringsAsFactors = F)

foobar3 <- tidyr::gather(foobar2, 
                         key="EvidenceClass",
                         value="Evidence",
                         c('Genetics',
                                 'Networks',
                                 'Expression',
                                 'ADPhenotype',
                                 'Metabolomic',
                                 'Proteomic',
                                 'ModelSystem'))

g <- ggplot2::ggplot(foobar3, 
                     ggplot2::aes(y=Gene,
                                  x=EvidenceClass,
                                  size=Evidence))
g <- g + ggplot2::geom_count()
g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
g


foobar4 <- foobar2
rownames(foobar4) <- foobar2$Gene
foobar4 <- foobar4[,-1]

pheatmap::pheatmap(foobar4)
