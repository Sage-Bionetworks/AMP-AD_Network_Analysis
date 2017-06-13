GenGraphViz <- function(x, pattern, Types,
                        sizeMin = 3, sizeMax = 20){
  
  #identifying subset of data to use
  print('Obtaining sub-graph based on brain region')
  In1 <- grep(pattern,Dat$ModuleNameFull1)
  In2 <- grep(pattern,Dat$ModuleNameFull2)
  In_int <- intersect(In1,In2)
  Dat <- Dat[In_int,]
  
  #annotating nodes by algorithm type
  print('Annotating nodes by algorithm type')
  ModNames <- c()
  ModType <- c()
  TypeName <- c()
  
  l <- length(Types)
  
  for (i in 1:l){
    temp <- union(grep(Types[i],Dat$ModuleNameFull1, value= TRUE),
                  grep(Types[i],Dat$ModuleNameFull2, value= TRUE))
    ModNames <- c(ModNames,temp)
    ModType <- c(ModType,rep(i,length(temp)))
    cat('No. of modules in',Types[i],length(temp),'\n')
    TypeName <- c(TypeName, rep(Types[i],length(temp)))
  }
  
  #creating visualization 
  print('Creating visualization')
  df <- data.frame(Name = ModNames, ModType = ModType, 
                   TypeName = TypeName)
  library('igraph')
  df2 <- data.frame(from = Dat$ModuleNameFull1, 
                    to = Dat$ModuleNameFull2,
                    weight = 1/(1e-50 + c(Dat$fisherOR)^1))
  net <- graph.data.frame(df2, vertices=df, directed=F) 
  
  #coloring nodes on the basis of module type 
  cmap <- rainbow(l, alpha=1) 
  V(net)$color <- cmap[V(net)$ModType]
  deg <- 1/(1e-50 + strength(net))
  V(net)$size <- sizeMin + deg/max(deg)*sizeMax
  l <- layout_with_fr(net, niter = 5000, grid = 'nogrid')
  plot(net, edge.arrow.size=.4,vertex.label=NA,  
       layout=l)
  legend(x=-1.5, y=-1.1,Types,pch=21,col="#777777", pt.bg=cmap,
         pt.cex=2, cex=.8, bty="n", ncol=1)
  
  #returning igraph object and layout
  RetDat <- list()
  RetDat$net <- net 
  RetDat$l <- l 
  return(RetDat)
  
  
}



GenClusteredViz <- function(Net,l){
  
  #loading libraries 
  library(igraph)
  library(MCL)
  
  #identifying clusters 
  print('Identifying clusters')
  Adj <- as_adj(Net, sparse = TRUE)
  Clusts <- mcl(Adj, allow1 = TRUE, addLoops = FALSE)
  cmap <- rainbow(Clusts$K, alpha=1)
  V(Net)$color <- cmap[Clusts$Cluster]
  
  cat('No. of clusters estimated',Clusts$K,'\n')
  
  for( i in 1:Clusts$K){
    cat('No. of modules in cluster',i,'is =',
        sum(Clusts$Cluster==i),'\n')
  }
    
  
  #generating the plot 
  print('Generating the plot')
  #plot(Net, edge.arrow.size=.4,vertex.label=NA,
   #    layout=layout_with_fr(Net, niter = 5000, grid = 'nogrid'))
  plot(Net, edge.arrow.size=.4,vertex.label=NA, layout = l)
  
  #adding legend 
  LegendNames <- c(1:Clusts$K)
  legend(x=-1.5, y=-1.1,LegendNames,pch=21,col="#777777", pt.bg=cmap,
         pt.cex=2, cex=.8, bty="n", ncol=6)
  
  ClustResults <- list()
  ClustResults$Net <- Net
  ClustResults$Clusts <- Clusts
  ClustResults$Adj <- Adj
  
  return(ClustResults)
  
}