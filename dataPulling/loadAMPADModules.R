synapseClient::synapseLogin()
mods <- list()

mods$DLPFC <- synapseClient::synTableQuery('SELECT * FROM syn9705614')@values
mods$IFG <- synapseClient::synTableQuery('SELECT * FROM syn9730668')@values
mods$STG <- synapseClient::synTableQuery('SELECT * FROM syn9730669')@values
mods$TCX <- synapseClient::synTableQuery('SELECT * FROM syn9730674')@values
mods$CER <- synapseClient::synTableQuery('SELECT * FROM syn9730675')@values
mods$PHG <- synapseClient::synTableQuery('SELECT * FROM syn9730672')@values
mods$FP <- synapseClient::synTableQuery('SELECT * FROM syn9737595')@values

#combine into a single df
addBrainRegionFxn <- function(x,y){
  x$brainRegion <- rep(y,nrow(x))
  return(x)
}

mods2 <- mapply(addBrainRegionFxn,mods,names(mods),SIMPLIFY=F)
allMods <- do.call(rbind,mods2)
allMods$ModuleNameFull <- paste0(allMods$ModuleName,allMods$brainRegion)
