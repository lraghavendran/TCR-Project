### This script takes Hamming's distance matrix for the samples and creates a network followed by cluster
### analysis on the edges.

### Install libraries
install.packages("dplyr", dependencies = TRUE)
install.packages("igraph", dependencies = TRUE)
install.packages("RColorBrewer", dependencies = TRUE)
install.packages("Polychrome", dependencies = TRUE)

### Load libraries
library('dplyr')
library('igraph')
library('RColorBrewer')
library('Polychrome')

### Define the input and output folders
data_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/data/"
output_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/output/"

### Load the Hamming's distance matrix for all the samples
input_file_name = paste0(data_folder,"ALL_Alignment_distance_matrix.csv")
ALL_distance <- read.csv(input_file_name, header = TRUE, row.names=1)

### Convert the distance dataframe to matrix of distances 
mat <- as.matrix(ALL_distance)

### Convert matrix of distances to matrix of similarities
mat <- max(mat)-mat

### convert everything to a scale of 0 to 1
mat <- mat/max(mat)

### shift the scale to [0,1] to convert negative correlatiosn to positive correlations
### not necessary for this dataset but will be helpful if we have a distance matrix with positive and negative values.
adj_mat <- (mat + 1)/2

### Emphirical filtering to get rid of low similarity correlations in this case below 0.6
### ((0.6+1)/2)**3 = 0.512
adj_mat <- adj_mat^3
adj_mat[adj_mat < 0.5] <- 0

### Remove self-loop maximum correlation value.
diag(adj_mat) <- 0

### create graph from the matrix using igraph functions
g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE)

### plot using a force-directed layout (fruchterman-reingold)
### Set seed for the network initialization down the code
set.seed(1)
coords_fr = layout.fruchterman.reingold(g, weights=E(g)$weight)

### To relate the nodes/vertices with the sample names
rownames_all <- row.names(ALL_distance)
nodes <- data.frame(id = 1:length(rownames_all), sample = rownames_all)

### Define groups to be able to name the nodes in the network results
nodes$group <- ifelse(grepl("M100[A-Z]",nodes$sample),'M100',
                      ifelse(grepl("M10[A-Z]",nodes$sample),'M10',
                             ifelse(grepl("M1[A-Z]",nodes$sample),'M1',
                                    'P')))
### Define group types to be able to color the nodes in the network results
nodes$group.type <- ifelse(grepl("M100[A-Z]",nodes$sample),'4',
                           ifelse(grepl("M10[A-Z]",nodes$sample),'3',
                                  ifelse(grepl("M1[A-Z]",nodes$sample),'2',
                                         '1')))
### Color types
### Control samples (P) - 1 - Yellow
### Treatment samples (M1) - 2 - Tomato (Red)
### Treatment samples (M10) - 3 - Green 
### Treatment samples (M100) - 4 - Blue

### Define the group.type to Color by sample type 
V(g)$group.type <- nodes$group.type
colrs <- c("yellow", "tomato", "green", "blue")
V(g)$color <- colrs[as.numeric(V(g)$group.type)]

### Visualize the network analysis plot. The randomness of the position initialization would mean the 
### clusters of samples would change their relative position, while the clusters remain intact.
plot(g, layout=coords_fr, edge.arrow.size=0, vertex.size=7, vertex.color=V(g)$color,
     vertex.label=nodes$sample, vertex.label.cex=.5, vertex.label.color="black")

### Cluster edge identification requires the weights to be greater than 0
E(g)$weight[E(g)$weight == 0] <- 0.1

### Identify clusters from the network analysis
ceb <- cluster_edge_betweenness(g)

### Modularity scores represents the connectiveness of all the edges within a community
modularity(ceb)

### How many clusters are identified from the cluster analysis
levels(factor(ceb$membership))

### Total number of clusters identified
tot_clus <- max(ceb$membership)

### Color the network analysis by cluster type
### Hexadecimal color specification for larger color palette
library(Polychrome)
Glasbey = glasbey.colors(tot_clus)
swatch(Glasbey)

### Define the memb.type to Color by cluster type
V(g)$memb.type <- ceb$membership

### Define colors
V(g)$color <- Glasbey[V(g)$memb.type]

### Visualize the network analysis plot. The randomness of the position initialization would mean the 
### clusters of samples would change their relative position, while the clusters remain intact.
### color by clusters

plot_file_name = paste0(output_folder,"network_plot_thresh_06.pdf")

# Step 1: Call the pdf command to start the plot
pdf(file = plot_file_name,   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
plot(g, layout=coords_fr, edge.arrow.size=0, vertex.size=7, vertex.color = adjustcolor(V(g)$color, alpha.f = .7),
     vertex.label=nodes$label, vertex.label.cex=.5, vertex.label.color="black")

#plot(ceb, g, layout=coords_fr, edge.arrow.size=0, vertex.size=7, vertex.color = adjustcolor(V(g)$color, alpha.f = .7),
     #vertex.label=nodes$label, vertex.label.cex=.5, vertex.label.color="black")

legend(-1.5, 1, levels(as.factor(V(g)$memb.type)) ,pch=21,
       col="#777777", pt.bg=adjustcolor(Glasbey, alpha.f = 0.7), pt.cex=2, cex=1, bty="n", ncol=2, x.intersp = 0.5)

# Step 3: Run dev.off() to create the file!
dev.off()

### collect sample names and cluster information (membership) of clusters in to separate dataframe
ceb_df <- data.frame(name = ceb$names,
                     membership = ceb$membership)

### Print out the cluster number for each sample in the data folder as we need this file for further analysis
out_file_name = paste0(data_folder,"sample_clus_membership.txt")

write.table(ceb_df, file = out_file_name, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, fileEncoding = "")
