#login to synapse
synapseClient::synapseLogin()

#download individual module definitions
individualModules <- synapseClient::synTableQuery("select * from syn10309369")@values

#download consensus module definitions
consensusModules <- synapseClient::synTableQuery("select * from syn10337531")@values

#combine
allModules <- rbind(individualModules,consensusModules)

#all modules
foo <- rSynapseUtilities::makeTable(allModules,
                                    "all module manifest August 16 2017",
                                    "syn2370594")
