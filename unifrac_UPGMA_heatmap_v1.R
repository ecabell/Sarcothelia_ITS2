# Combined Phylogenetic Analysis - Version 1 (Panel A above, Panel B centered below)
# R script to recreate the unifrac distance analysis with dendrogram and heatmap

# Load required libraries
library(ggplot2)
library(pheatmap)
library(cluster)
library(dendextend)
library(gridExtra)
library(grid)
library(RColorBrewer)

# Read the unifrac distance file
# Note: Adjust the file path as needed
unifrac_file <- "~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist"

# Read the data
data_lines <- readLines(unifrac_file)

# Parse the data (skip the ID column - column 2)
sample_names <- c()
distance_data <- list()

for (i in 1:length(data_lines)) {
  line <- data_lines[i]
  if (nchar(line) > 0) {
    parts <- strsplit(line, "\t")[[1]]
    sample_name <- parts[1]
    # Skip parts[2] which is the ID column
    distances <- as.numeric(parts[3:length(parts)])
    
    sample_names <- c(sample_names, sample_name)
    distance_data[[i]] <- distances
  }
}

# Create distance matrix
distance_matrix <- do.call(rbind, distance_data)
rownames(distance_matrix) <- sample_names
colnames(distance_matrix) <- sample_names

# Make matrix symmetric and ensure diagonal is zero
distance_matrix_sym <- (distance_matrix + t(distance_matrix)) / 2
diag(distance_matrix_sym) <- 0

# Convert to distance object for clustering
dist_obj <- as.dist(distance_matrix_sym)

# Perform hierarchical clustering
hclust_result <- hclust(dist_obj, method = "average")

# Create dendrogram
dend <- as.dendrogram(hclust_result)

# Set up the plotting layout
png("unifrac_phylogenetic_analysis_v1_R.png", width = 16, height = 12, units = "in", res = 300)

# Create layout matrix for custom positioning
layout_matrix <- matrix(c(1, 1, 1,
                         0, 0, 0,
                         0, 2, 0), nrow = 3, byrow = TRUE)
layout(layout_matrix, heights = c(1, 0.1, 1.2))

# Set margins
par(mar = c(5, 4, 4, 2))

# Panel A: Dendrogram (top)
plot(dend, 
     main = "Neighbor-Joining Phylogenetic Tree of ITS2 Profiles Weighted UniFrac Distances",
     xlab = "ITS2 Profiles",
     ylab = "UniFrac Distance",
     cex.main = 1.4,
     cex.lab = 1.2,
     cex.axis = 1.0)
grid(ny = NULL, nx = NA, col = "lightgray", lty = "dotted")

# Add panel label A
mtext("A", side = 3, line = 2, at = par("usr")[1] - 0.1 * diff(par("usr")[1:2]), 
      cex = 1.6, font = 2)

# Panel B: Heatmap (bottom center)
par(mar = c(8, 8, 4, 4))

# Get the order from clustering
cluster_order <- order.dendrogram(dend)
ordered_matrix <- distance_matrix_sym[cluster_order, cluster_order]
ordered_names <- sample_names[cluster_order]

# Create color palette
colors <- colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100)

# Create heatmap
image(1:nrow(ordered_matrix), 1:ncol(ordered_matrix), ordered_matrix,
      col = colors,
      xlab = "", ylab = "",
      main = "UniFrac Distance Matrix Heatmap\n(Ordered by Phylogenetic Relationships)",
      cex.main = 1.4,
      axes = FALSE)

# Add axis labels
axis(1, at = 1:length(ordered_names), labels = ordered_names, 
     las = 2, cex.axis = 0.9)
axis(2, at = 1:length(ordered_names), labels = ordered_names, 
     las = 2, cex.axis = 0.9)

# Add text annotations
for (i in 1:nrow(ordered_matrix)) {
  for (j in 1:ncol(ordered_matrix)) {
    text_color <- if (ordered_matrix[i, j] > max(ordered_matrix) * 0.6) "white" else "black"
    text(j, i, sprintf("%.3f", ordered_matrix[i, j]), 
         cex = 0.8, col = text_color)
  }
}

# Add colorbar
colorbar_x <- par("usr")[2] + 0.05 * diff(par("usr")[1:2])
colorbar_y <- seq(par("usr")[3], par("usr")[4], length.out = 100)
for (i in 1:99) {
  rect(colorbar_x, colorbar_y[i], colorbar_x + 0.02 * diff(par("usr")[1:2]), 
       colorbar_y[i+1], col = colors[i], border = NA)
}

# Add colorbar labels
colorbar_labels <- seq(0, max(distance_matrix_sym), length.out = 6)
colorbar_positions <- seq(par("usr")[3], par("usr")[4], length.out = 6)
text(colorbar_x + 0.04 * diff(par("usr")[1:2]), colorbar_positions, 
     sprintf("%.3f", colorbar_labels), cex = 0.8)
text(colorbar_x + 0.02 * diff(par("usr")[1:2]), par("usr")[4] + 0.05 * diff(par("usr")[3:4]), 
     "UniFrac Distance", cex = 1.0, srt = 90)

# Add panel label B
mtext("B", side = 3, line = 2, at = par("usr")[1] - 0.1 * diff(par("usr")[1:2]), 
      cex = 1.6, font = 2)

dev.off()

print("Version 1 R plot saved as unifrac_phylogenetic_analysis_v1_R.png")
