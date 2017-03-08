winsorizeAndScale <- function(x){
  library(dplyr)
  library(utilityFunctions)
  x <- t(x)
  x <- apply(x,2,utilityFunctions::winsorize)
  x <- scale(x)
  return(x)
}