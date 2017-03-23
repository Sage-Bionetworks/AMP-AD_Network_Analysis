synapseClient::synapseLogin()
foo <- synapseClient::synGet('syn8268669')
load(foo@filePath)
rosmap <- igraph::graph_from_adjacency_matrix(bicNetworks$network,
                                              mode='undirected')


collateGraphStatistics <- function(graph){
  model <- list()
  cat('Computing Degree...\n')
  model$degree <- igraph::degree(graph)
  #model$alpha_centrality <- igraph::alpha_centrality(graph,tol=1e-14)
  cat('Computing Authority Score...\n')
  model$authority_score <- igraph::authority_score(graph)$vector
  
  cat('Computing Closeness...\n')
  model$closeness <- igraph::closeness(graph)
  
  cat('Computing Eccentricity...\n')
  model$eccentricity <- igraph::eccentricity(graph)
  
  cat('Computing Eigenvector Centrality...\n')
  model$eigen_centrality <- igraph::eigen_centrality(graph)$vector
  
  cat('Computing Betweeness Centrality...\n')
  model$centr_betw <- igraph::betweenness(graph)
  
  cat('Computing Transitivity...\n')
  model$transitivity <- igraph::transitivity(graph,
                                             type='undirected')
  
  return(model)
}

rosmapNetworkStats <- collateGraphStatistics(rosmap)
rosmapNetworkStats <- data.frame(rosmapNetworkStats,
                                 stringsAsFactors=F)

rosmapNetworkStats <- dplyr::filter(rosmapNetworkStats,degree>0)

str(rosmapNetworkStats)
hist(rosmapNetworkStats$degree)
hist(rosmapNetworkStats$authority_score$vector)
n1 <- names(sort(rosmapNetworkStats$authority_score$vector,decreasing=T)[1:5])
n1 <- names(sort(rosmapNetworkStats$degree,decreasing=T)[1:5])

utilityFunctions::convertEnsemblToHgnc(n1)
