# Run Consensus module notebook

First, query the Network manifest fileview for all of the relevant modules for a given brain region (e.g. for frontal pole)

`moduleManifest <- synapseClient::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) AND ( analysisType = 'moduleIdentification' ) AND ( method = 'CFinder' OR method = 'fast_greedy' OR method = 'GANXiS' OR method = 'infomap' OR method = 'label_prop' OR method = 'linkcommunities' OR method = 'linkcommunities' OR method = 'louvain' OR method = 'spinglass' OR method = 'walktrap' ) AND (tissueOfOrigin = 'frontalPole') )")`

next, identify the parentid for these files:
`names(table(moduleManifest$values$parentId))`



Next we use this function to build the submission scripts:

https://github.com/th1vairam/metanetworkSynapse/blob/e6237e99efaf46b599dac825d4808638810b9f09/makeConsensusModuleSubmission.R





next, we use these functions: https://github.com/th1vairam/metanetwork/blob/9c63003cd466c52816dc79e7757335fb27175f24/R/findModules.consensusKmeans.R

and https://github.com/th1vairam/metanetworkSynapse/blob/e6237e99efaf46b599dac825d4808638810b9f09/buildConsensusModules.R

to test the consensus module building.





qsub -cwd -V -pe mpi 4 -o ./submission.scripts/syn8340017.syn8340019.syn8379809.kmeans.out -e ./submission.scripts/syn8340017.syn8340019.syn8379809.kmeans.err ./submission.scripts/syn8340017.syn8340019.syn8379809.kmeans.sh