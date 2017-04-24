synapseClient::synapseLogin()

#see status of bic networks

library(dplyr)
bicNetworkManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( method = 'bic' ) AND ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) )")
View(bicNetworkManifest@values)

consensusModuleManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'consensusModuleIdentification' ) )")
View(consensusModuleManifest@values)

moduleManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'moduleIdentification' ) AND ( method = 'CFinder' OR method = 'fast_greedy' OR method = 'GANXiS' OR method = 'infomap' OR method = 'label_prop' OR method = 'linkcommunities' OR method = 'linkcommunities' OR method = 'louvain' OR method = 'spinglass' OR method = 'walktrap' ) )")
View(moduleManifest@values)
table(moduleManifest@values$tissueOfOrigin)


consortiaModules <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'consensusModuleIdentification' OR analysisType = 'moduleIdentification')  AND (method = 'kmeans' OR method = 'wina' OR method = 'speakeasy' OR method = 'megena' OR method = 'wgcna')) AND (columnScaled = 'TRUE')")
View(consortiaModules@values)
