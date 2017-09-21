#this script reads the information from a csv file 
# and creates graph using it 

setwd('Documents/SageStuff/')

#subset selection 
Dat <- read.csv('Job-38889986948603617091033777.csv')
Dat <- data.frame(Dat)
pattern <- "DLPFC"
In1 <- grep(pattern,Dat$ModuleNameFull1)
In2 <- grep(pattern,Dat$ModuleNameFull2)
In_int <- intersect(In1,In2)
Dat <- Dat[In_int,]

#obtaining the type for each node 
Types <- c('consen','megen','metan','speakE','wina')
ModNames <- c()
ModType <- c()
TypeName <- c()

for (i in 1:5){
  temp <- union(grep(Types[i],Dat$ModuleNameFull1, value= TRUE),grep(Types[i],Dat$ModuleNameFull2, value= TRUE))
  ModNames <- c(ModNames,temp)
  ModType <- c(ModType,rep(i,length(temp)))
  print(length(temp))
  TypeName <- c(TypeName, rep(Types[i],length(temp)))
}

df <- data.frame(Name = ModNames, ModType = ModType, TypeName = TypeName)

library('igraph')
df2 <- data.frame(from = Dat$ModuleNameFull1, to = Dat$ModuleNameFull2, weight = c(Dat$fisherOR))
net <- graph.data.frame(df2, vertices=df, directed=F) 

#coloring nodes on the basis of module type 
cmap <- rainbow(5, alpha=1) 
V(net)$color <- cmap[V(net)$ModType]
deg <- strength(net)
V(net)$size <- 3 + deg/max(deg)*20
plot(net, edge.arrow.size=.4,vertex.label=NA,  layout=layout_with_fr(net, niter = 10000))
legend(x=-1.5, y=-1.1,c("Consensus","Megena", "Meta","SpeakEasy","Wina"),pch=21,col="#777777", pt.bg=cmap, pt.cex=2, cex=.8, bty="n", ncol=1)
