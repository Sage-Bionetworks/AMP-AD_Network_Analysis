####ROSMAP
synapseClient::synapseLogin()
#get manifest
foo <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( analysisType = 'moduleIdentification' ) AND ( method = 'GANXiS' OR method = 'fast_greedy' OR method = 'CFinder' OR method = 'linkcommunities' OR method = 'label_prop' OR method = 'louvain' OR method = 'infomap' OR method = 'spinglass' OR method = 'walktrap' ) AND ( tissueTypeAbrv = 'DLPFC' ) )")@values

#pull data
bar <- lapply(foo$id,synapseClient::synGet)
modules <- lapply(bar,function(x){
  library(dplyr)
  synapseClient::getFileLocation(x) %>%
    data.table::fread(data.table=F)})

baz <- sapply(modules,function(x) x$moduleLabel)
baz[1:5,1:5]
rownames(baz) <- modules[[1]]$Gene.ID
colnames(baz) <- foo$method
baz <- data.frame(baz)
#str(baz)
#baz2 <- model.matrix(~label_prop + fast_greedy + walktrap + spinglass + infomap + louvain + GANXiS + linkcommunities + CFinder-1,baz)

baz2 <- model.matrix( ~ . -1, data=baz, contrasts.arg = 
                lapply(data.frame(baz[,sapply(data.frame(baz), is.factor)]),
                       contrasts, contrasts = FALSE))

baz3 <- metanetwork::findModules.consensusKmeans(baz2[,sample(1:ncol(baz2))])

library(data.table)
partition.adj = lapply(foo$id, function(id){
  fread(synGet(id)@filePath, header = T, data.table =F)
})
names(partition.adj) = paste0('Method', 1:length(partition.adj))

partition.adj = mapply(function(mod, method){
  mod = mod %>%
    dplyr::select(Gene.ID, moduleNumber) %>%
    dplyr::mutate(value = 1,
                  moduleNumber = paste0(method,moduleNumber)) %>%
    tidyr::spread(moduleNumber, value)
}, partition.adj, names(partition.adj), SIMPLIFY = F) %>%
  join_all()
partition.adj[is.na(partition.adj)] = 0
rownames(partition.adj) = partition.adj$Gene.ID
partition.adj$Gene.ID = NULL

#### Compute consensus modules using specified algorithm ####
# Get a specific algorithm
findModules.algo = switch (cons.method,
                           kmeans = metanetwork::findModules.consensusKmeans)

# Compute consensus modules
mod <- findModules.algo(data.matrix(partition.adj), min.module.size = 20, usepam = FALSE)
mod <- metanetwork::findModules.consensusKmeans(data.matrix(partition.adj), min.module.size = 20, usepam= FALSE)

set.seed(1)
baz3 = fpc::pamk(baz2, krange = 2:30, usepam = FALSE, rngR = TRUE)
table(baz3$pamobject$cluster)
set.seed(1)
baz3 = fpc::pamk(data.matrix(partition.adj), krange = 2:30, usepam = FALSE, rngR = TRUE)
table(baz3$pamobject$cluster)


set.seed(1)
baz3 = fpc::pamk(baz2, krange = 2:30, usepam = FALSE, rngR = TRUE)
table(baz3$pamobject$cluster)
set.seed(1)
baz3 = fpc::pamk(baz2, krange = 2:30, usepam = FALSE, rngR = TRUE,criterion = "multiasw")
table(baz3$pamobject$cluster)
set.seed(1)
baz3 = fpc::pamk(baz2[sample(1:nrow(baz2),nrow(baz2)),], krange = 2:30, usepam = FALSE, rngR = TRUE)
table(baz3$pamobject$cluster)

set.seed(1)
baz3 = fpc::pamk(baz2[sample(1:nrow(baz2),nrow(baz2)),], krange = 2:30, usepam = FALSE, rngR = TRUE,criterion = "multiasw")
table(baz3$pamobject$cluster)


set.seed(12345);clust1 = fpc::pamk(baz2, usepam = F, k = 3:10, seed = 12345, criterion = 'ch', rngR = T, sampsize = 50);table(clust1$pamobject$cluster)
set.seed(12345);clust1 = fpc::pamk(baz2[sample(1:nrow(baz2),nrow(baz2)),], usepam = F, k = 3:10, seed = 12345, criterion = 'ch', rngR = T, sampsize = 50);table(clust1$pamobject$cluster)

baz5 <- baz2[sample(1:nrow(baz2),1000),]

set.seed(12345);clust1 = fpc::pamk(baz5, usepam = T, k = 3:10);sort(table(clust1$pamobject$cluster))
set.seed(12345);clust1 = fpc::pamk(baz5[sample(1:nrow(baz5),1000),], usepam = T, k = 3:10);sort(table(clust1$pamobject$cluster))

modDf <- sapply(modules,function(x) x$moduleNumber)
rownames(modDf) <- modules[[1]]$Gene.ID
colnames(modDf) <- foo$method
#####try out clue consensus
getClue <- function(x,y){
  names(x) <- y
  return(clue::as.cl_partition(x))
}
rosmap_modules <- apply(modDf,2,getClue,rownames(modDf))
he_rosmap <- clue::cl_consensus(rosmap_modules,method = 'HE')
dwh_rosmap <- clue::cl_consensus(rosmap_modules,method = 'DWH')

dwh_rosmap_hard <- clue::as.cl_hard_partition(dwh_rosmap)
table(se_rosmap_hard$.Data)
sum(clue::cl_dissimilarity(rosmap_modules, se_rosmap, "comem") ^ 2)

##### new strategy
#run simple kmeans (try n = 1, 10, 100, 1000) for a given k
#use as initialization for flexmix
#run flexmix
#find minimum aic over path
#do deep dive on minimum
#return cluster definitions

####flex mix
install.packages("flexmix")
set.seed(1)
k1 <- 5
kmeRes <- stats::kmeans(baz2,k1,nstart = 1)
table(kmeRes$cluster)

multivariateNormalDensity <- function(x,mu,sigma){
  library(dplyr)
  dn<-0.5*log(2*pi*det(sigma)) - 0.5*t(x-mu)%*%solve(sigma)%*%(x-mu) %>% as.numeric
  return(dn)
}



computeKmeansScores <- function(kmat,kvec){
  model <- list()
  muVec <- apply(kmat,2,mean)
  sdVec <- apply(kmat,2,sd)
  varMat <- diag(sdVec^2)
  dVec <- apply(kmat,1,multivariateNormalDensity,muVec,varMat)
  #countVec <- table(kvec)
  #probVec <- table(kmeRes$cluster)/length(kmeRes$cluster)
  #countVec <- table(kmeRes$cluster)
  model$loglik <- sum(dVec)
  
  model$bic <- -2*(model$loglik) + ncol(kmat)*log(length(kvec))
  model$aic <- -2*(model$loglik) + 2*ncol(kmat)
  return(model)
}
k1 <- 2
kmeRes <- stats::kmeans(baz2,k1,nstart = 1)
kmeansDf <- data.frame(kmeRes$cluster)
kmeansDf$kmeRes.cluster <- as.factor(kmeansDf$kmeRes.cluster)
kmat <- model.matrix( ~ . -1, data=kmeansDf, contrasts.arg = 
                        lapply(data.frame(kmeansDf[,sapply(data.frame(kmeansDf), is.factor)]),
                               contrasts, contrasts = FALSE))
computeKmeansScores(kmat,kmeRes$cluster)



install.packages('BayesLCA')
kmeRes <- stats::kmeans(baz2,2,nstart = 1)

res2 <- BayesLCA::blca(t(baz2[1:5000,]),3,method='vb',start.vals = "across",restarts=1,iter = 1e3, conv=1e-9)
res3 <- BayesLCA::blca(t(baz2),3,method='em',restarts=1)

score20 <- computeKmeansScores(kmeRes$cluster)
score20

kmat <- model.matrix( ~ . -1, data=kmeansDf, contrasts.arg = 
                lapply(data.frame(kmeansDf[,sapply(data.frame(kmeansDf), is.factor)]),
                       contrasts, contrasts = FALSE))
m1 <- flexmix::flexmix(baz2 ~ 1,
                       k = 11,
                       model = flexmix::FLXMCmvbinary())
a5 <- apply(m1@posterior$scaled,
            1,
            which.max)
table(a5,kmeRes$cluster)
#m1 <- flexmix::initFlexmix(baz2~1, k = 2:30, model = FLXMCmvbinary(),nrep = 1)
#m1 <- flexmix::flexmix(baz2 ~ 1, k = 4,model = flexmix::FLXMCmvbinary())
AIC(m1)
a5 <- m1@cluster
#a5<-apply(m1@models[[3]]@posterior$scaled,1,which.max)
#table(a5)
#image(m1@posterior$scaled)

data("Nclus")
m1 <- flexmix(Nclus ~ 1, k = 4, model = mymclust())
R> summary(m1)


library(LEA)
#apply(baz2,1,function(x,fileName) {cat(x,'\n',sep='',file=fileName,append=T)},'test.geno')
#cat(baz2,sp)
#foobar <- LEA::pca('test.geno')
baz2save <- baz2
#baz2 <- matrix(rbinom(1e5,2,.5),1e3,1e2)
baz2 <- baz2save
svd3 <- svd((apply(baz2,2,scale)))
#install.packages('RMTstat')
library(RMTstat)
nprime <- function(numberOfCluster,eigenvals){
  a1 <-((numberOfCluster+1)*((sum(eigenvals))^2))
  a2 <- a1/((numberOfCluster-1)*((sum(eigenvals^2))) -((sum(eigenvals))^2))
  model <- list()
  model$nprime <- a2
  model$sigma2 <- sum(eigenvals)/((numberOfCluster-1)*a2)
  return(model)
}
npr<-nprime(ncol(baz2),svd3$d^2)
computeTWstat <- function(i,baz2,svd3){
  mprime <- ncol(baz2)-1-i
  nprime <- nrow(baz2)-1-i
  n <- nrow(baz2)-i
  m <- ncol(baz2)-i
  ell <- (mprime*(svd3$d^2)[i+1])/(sum((svd3$d[(i+1):mprime]^2)))
  mu <- ((sqrt(n-1) + sqrt(m))^2)/(n)
  sig <- ((sqrt(n-1)+sqrt(m))/(n))*((1/(sqrt(n-1))+1/sqrt(m))^(1/3))
  #lambda1 <- (svd3$d^2)[i+1]
  x <- (ell-mu)/sig
  return(x)
}
twstats <- sapply(1:270,computeTWstat,(baz2),svd3)
#RMTstat::dtw(x,beta=4,log=F)

kmeansAIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

computeAICfxn <- function(k,baz2,nstartPar){
  metanetwork2 <- kmeans(baz2,k,nstart=nstartPar);
  foob<-(kmeansAIC(metanetwork2))
  cat('k:',k,'aic',foob,'\n')
  return(foob)
}
set.seed(1)
aics <- sapply(2:60,computeAICfxn,baz2,1)
aicsRefined <- sapply(10:35,computeAICfxn,baz2,10)
aicsRefined2 <- sapply(17:28,computeAICfxn,baz2,100)
plot(2:60,aics,'l',ylim=c(5e4,1e5),xlim=c(1,60),col='blue',lwd=2)
par(new=T)
plot(10:35,aicsRefined,'l',ylim=c(5e4,1e5),xlim=c(1,60),col='purple',lwd=2)
par(new=T)
plot(17:28,aicsRefined2,'l',ylim=c(5e4,1e5),xlim=c(1,60),col='red',lwd=2)
plot(aics)

set.seed(1)
t1 <- kmeans(baz2,24,nstart=1000)
set.seed(2)
t2 <- kmeans(baz2,24,nstart=1000)
foo1 <- table(t1$cluster,t2$cluster)
computePercentageDiscordant <- function(x){
  ntot <- sum(x)
  x <- apply(x,1,function(y){a <- y;a[which.max(y)]<-0;return(a)})
  return(sum(x)/ntot)
}
computePercentageDiscordant(foo1)

mean(rowSums(foo1>0))

plot(rowSums(foo1>0))

####NEW STRATEGY: use AIC to calculate the number of clusters supported.


set.seed(1)
metanetwork172 <- kmeans(baz2,172,nstart = 100)
sortedCounts <- metanetwork172$cluster %>%
  table %>%
  sort
pairs(svd3$u[,1:4],col=rainbow(172)[metanetwork172$cluster])
