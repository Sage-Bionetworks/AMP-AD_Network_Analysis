#network comparison script
library(synapseClient)
synapseLogin()
foo <- synQuery('select * from file where method==\'bic\' and projectId==\'syn2370594\'')

bar <- dplyr::select(foo,file.name,file.id,file.versionComment)

