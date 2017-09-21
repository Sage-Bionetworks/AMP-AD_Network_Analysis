l <- list()
source('TestConsensus_Sumit.R')
l$CBE <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/CBE_CCPivot_Consensus.rds'),5)
l$DLPFC <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/DLPFC_CCPivot_Consensus.rds'),5)
l$FP <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/FP_CCPivot_Consensus.rds'),5)
l$IFG <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/IFG_CCPivot_Consensus.rds'),5)
l$PHG <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/PHG_CCPivot_Consensus.rds'),5)
l$STG <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/STG_CCPivot_Consensus.rds'),5)
l$TCX <- ClusterLowFreqConsensus(
  readRDS('ConsensusMods/CC-Pivot/TCX_CCPivot_Consensus.rds'),5)

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'CBE',"' ) )")
CBE_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'DLPFC',"' ) )")
DLPFC_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'FP',"' ) )")
FP_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'IFG',"' ) )")
IFG_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'PHG',"' ) )")
PHG_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'STG',"' ) )")
STG_Mods <- synapseClient::synTableQuery(queryString)@values

queryString <- paste0("SELECT * FROM ",manifestId," WHERE ( ( brainRegion = '",'TCX',"' ) )")
TCX_Mods <- synapseClient::synTableQuery(queryString)@values


CBE <- GenConsMods(l$CBE, 'consensus', CBE_Mods , 'CBE')
DLPFC <- GenConsMods(l$DLPFC, 'consensus', DLPFC_Mods, 'DLPFC')
FP <- GenConsMods(l$FP, 'consensus', FP_Mods, 'FP')
IFG <- GenConsMods(l$IFG, 'consensus', IFG_Mods, 'IFG')
PHG <- GenConsMods(l$PHG, 'consensus', PHG_Mods, 'PHG')
STG <- GenConsMods(l$STG, 'consensus', STG_Mods, 'STG')
TCX <- GenConsMods(l$TCX, 'consensus', TCX_Mods, 'TCX')

AllMods <- rbind(CBE,DLPFC,FP,IFG,PHG,STG,TCX)
head(AllMods)
source('AMP-AD_Network_Analysis/UtilityFiles/SumitUtilFcns.R')
makeTable(AllMods,'Consensus modules all regions CombReps','syn5569099')
