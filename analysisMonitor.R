synapseClient::synapseLogin()

#see status of bic networks

library(dplyr)
bicNetworkManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( method = 'bic' ) AND ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) )")
View(bicNetworkManifest@values)

consensusModuleManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'consensusModuleIdentification' ) )")
View(consensusModuleManifest@values)

consortiaModules <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'consensusModuleIdentification' OR analysisType = 'moduleIdentification')  AND (method = 'kmeans' OR method = 'wina' OR method = 'speakeasy' OR method = 'megena' OR method = 'wgcna')) AND (columnScaled = 'TRUE')")
