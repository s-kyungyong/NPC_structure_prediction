library(igraph)

# Read data
data <- read.table(file = 'NPC_prediction_scores.txt', sep = '\t', header = TRUE)

# Generate a network
graph_data <- data[data$ipTM > 0.75 | data$interfaceScore > 0.5, ]
graph <-  graph.data.frame(graph_data, directed = FALSE)

# Color edges
edge_colors <- ifelse(graph_data$ipTM > 0.75 & graph_data$interfaceScore > 0.5, "green",
               ifelse(graph_data$ipTM > 0.75 , "grey",
               ifelse(graph_data$interfaceScore  > 0.5, "black", 
                                    NA)))
E(graph)$color <- edge_colors

# Plot the network graph with enhanced aesthetics
set.seed(321)
plot(graph,
     edge.width = 2,
     edge.arrow.size = 0.5,
     edge.color = E(graph)$color,
     vertex.color = "lightblue",
     vertex.size = 7,
     vertex.label = V(graph)$name,
     vertex.label.cex = 0.8,
     vertex.frame.color = "gray",
     vertex.frame.width = 2,
     layout = layout_with_fr(graph, niter = 1000, area = 4000),
     main = "Network Graph",
     sub = "Gene Connections",
     margin = 0.1)
