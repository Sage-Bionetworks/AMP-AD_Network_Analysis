synapseClient::synapseLogin()
allMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ","syn10146524"))@values
consensusMods <- synapseClient::synTableQuery(paste0("SELECT * FROM ","syn10158502"))@values
allMods <- rbind(allMods,consensusMods)
rSynapseUtilities::makeTable(allMods,"full module manifest July 6 2017",'syn2370594')
