synapseClient::synapseLogin()
rosmapClinicalObj <- synapseClient::synGet('syn3191087')
rosmapUncensoredAgesObj <- synapseClient::synGet('syn7116000')
rosmapIdMapObj <- synapseClient::synGet('syn3382527')
rosmapCogDecline1Obj <- synapseClient::synGet('syn6182375')
rosmapCogDecline2Obj <- synapseClient::synGet('syn6182376')
rosmapMetanetworkObj <- synapseClient::synGet('syn8268669')

####thanneer's code for fixing ids
#https://github.com/th1vairam/Brain_Reg_Net/blob/ad12b544d6d6c23be2bf94eed53a6b4f75d154d1/code/Rmd/ROSMAP_REPROCESSED.Rmd


