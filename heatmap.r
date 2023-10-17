library(pheatmap)
library(RColorBrewer)

# Read data
data <- read.table(file = 'NPC_prediction_scores.txt', sep = '\t', header = TRUE)

# Subset data
data_subset <- subset(data, ipTM > 0.75)

# Create a new matrix new row names
concat_names <- paste(data_subset$Pair1, data_subset$Pair2, sep = " ")
data_new <- matrix( data = data_subset$ipTM, ncol = 1)
rownames(data_new) <- concat_names

# Create heatmaps
cols <- c(colorRampPalette(brewer.pal(9, "YlOrRd"))(100))
pheatmap(data_new, color = cols,
         cluster_row  = FALSE, cluster_cols = FALSE )
