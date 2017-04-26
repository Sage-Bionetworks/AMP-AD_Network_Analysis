synapseClient::synapseLogin()

#get synIds for gene expression variables
geneExpressionDataManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( dataSubType = 'residualGeneExpForNetAnlz' ) AND ( normalizationType = 'CQN' ) )")

#get synIds for covariates
covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

