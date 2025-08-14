# Install and load necessary packages if you haven't already
library(ape)
library(pheatmap)

# Read the raw file without header
raw_data <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract the sample names
sample_names <- raw_data[[1]]

# Remove the first two columns (sample names and UIDs)
dist_data <- raw_data[, -(1:2)]

# Convert to matrix
dist_matrix <- as.matrix(dist_data)

# Assign row and column names
rownames(dist_matrix) <- sample_names
colnames(dist_matrix) <- sample_names

# Plot heatmap
pheatmap(dist_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         display_numbers = TRUE,         # Show distance values
         number_format = "%.3f",         # Optional: format to 3 decimal places
         main = "Symbiodiniaceae ITS2 Profiles Unifrac Distance Heatmap",
         fontsize = 12,
         cellwidth = 50,
         cellheight = 40)

# Read the aligned DNA sequences (FASTA format)
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# Compute pairwise genetic distances
genetic_dist <- dist.dna(alignment, model = "raw")  # Use "JC69", "K80", etc. for evolutionary models

# Convert to a matrix and inspect
genetic_matrix <- as.matrix(genetic_dist)
print(round(genetic_matrix, 4))

# Read mapping
map <- read.csv("~/Documents/Zoox PCR/meta.csv", stringsAsFactors = FALSE)

# Replace labels
new_labels <- setNames(map$ITS2_Type, map$UID)

# Apply relabeling to matrix
rownames(genetic_matrix) <- new_labels[rownames(genetic_matrix)]
colnames(genetic_matrix) <- new_labels[colnames(genetic_matrix)]

#heatmap
pheatmap(genetic_matrix,
         display_numbers = TRUE,
         number_format = "%.3f",
         main = "FASTA-Based Symbiodiniaceae ITS2 Type Genetic Distance Heatmap",
         fontsize = 12,
         cellwidth = 50,
         cellheight = 40)


# Compute corrected genetic distances
genetic_dist <- dist.dna(alignment, model = "JC69")  # Use "JC69", "K80", etc. for evolutionary models

# Convert to a matrix and inspect
genetic_matrix <- as.matrix(genetic_dist)
print(round(genetic_matrix, 4))

# Read mapping
map <- read.csv("~/Documents/Zoox PCR/meta.csv", stringsAsFactors = FALSE)

# Replace labels
new_labels <- setNames(map$ITS2_Type, map$UID)

# Apply relabeling to matrix
rownames(genetic_matrix) <- new_labels[rownames(genetic_matrix)]
colnames(genetic_matrix) <- new_labels[colnames(genetic_matrix)]

#heatmap
pheatmap(genetic_matrix,
         display_numbers = TRUE,
         number_format = "%.3f",
         main = "FASTA-Based Symbiodiniaceae ITS2 Type Genetic Distance Heatmap",
         fontsize = 12,
         cellwidth = 50,
         cellheight = 40)


