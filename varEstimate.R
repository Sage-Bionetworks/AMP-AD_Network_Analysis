###variance estimation
library(synapseClient)
synapseLogin()

foo <- synGet('syn7984097')

library(data.table)

rosmapExpression <- fread(foo@filePath,data.table=F)

barObj <- synGet('syn8018344')

rosmapClinical <- fread(barObj@filePath,data.table=F)

rownames(rosmapExpression) <- rosmapExpression[,1]
rosmapExpression <- rosmapExpression[,-1]
rosmapExpression <- t(rosmapExpression)
varCompWrapperFunction <- function(outcome,features){
  library(varComp)
  ####features have to be in observations as rows and features as column format
  kernel <- cor(t(features),
                use = 'pairwise.complete.obs')
  alternativeModel <- varComp::varComp(outcome~1,
                                       varcov=kernel)
  nullModel <- varComp::varComp(outcome~1)
  
  model <- list()
  model$totalVariance <- as.numeric(var(outcome))
  model$varianceExplained <- as.numeric(alternativeModel$varComps[1])
  model$percentVarianceExplained <- model$varianceExplained/model$totalVariance
  model$errorVarianceExplained <- as.numeric(alternativeModel$sigma2)
  model$LRT <- 2*(logLik(alternativeModel)[1]-logLik(nullModel)[1])
  
  pVarianceLRT <- function(x){
    x <- round(x,8)
    if(x==0){
      return(0.5)
    }else{
      return(0.5*pchisq(x,1,lower.tail=F))
    }
  }
  
  model$pvalue <-pVarianceLRT(model$LRT)
  return(model)
}
bar4<-fendR::varCompWrapperFunction(rosmapClinical$pmi,rosmapExpression)

