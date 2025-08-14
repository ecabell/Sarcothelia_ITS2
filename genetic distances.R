# Install required packages if not already installed
install.packages("ape")
install.packages("seqinr")

# Load libraries
library(ape)
library(seqinr)

# Load aligned sequences from FASTA file
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# Calculate pairwise distances using the p-distance (proportion of nucleotide differences)
distance_matrix <- dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)

# Convert to matrix for display
distance_matrix <- as.matrix(distance_matrix)

# View first few rows
print(round(distance_matrix, 4))

# Compute average pairwise genetic distance (excluding self-comparisons)
avg_distance <- mean(distance_matrix[upper.tri(distance_matrix)])
cat("Average pairwise genetic distance:", round(avg_distance, 4), "\n")


# Optional: Save to CSV
write.csv(distance_matrix, "genetic_distance_matrix.csv")

# Compute distance matrix
dist_matrix <- dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)

# Build a tree (e.g., Neighbor-Joining)
tree <- nj(dist_matrix)

# Plot tree
plot(tree, main = "Phylogenetic Tree with Genetic Distances")
nodelabels(round(tree$edge.length, 4), frame = "n", cex = 0.6)

library(gplots)
heatmap(as.matrix(dist_matrix), symm = TRUE, margins = c(10,10))

# Load aligned DNA sequences
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# Compute uncorrected p-distance
p_dist <- dist.dna(alignment, model = "raw")

# View as matrix
as.matrix(p_dist)

# Check names
rownames(alignment)

# Read mapping
map <- read.csv("~/Documents/Zoox PCR/meta.csv", stringsAsFactors = FALSE)

# Replace labels
new_labels <- setNames(map$ITS2_Type, map$UID)


# Subset to the two sequences of interest
pair <- alignment[c("22385", "760355"), ]  # replace with actual names

# Compute pairwise p-distance
dist <- dist.dna(pair, model = "raw")
as.matrix(dist)


