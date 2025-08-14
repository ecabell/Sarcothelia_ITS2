# Load required packages
library(vegan)
library(readr)

# Load distance matrix
dist_raw <- read_delim("~/Downloads/20250603T195144_cabeler1/between_sample_distances/A/20250605T195712_braycurtis_sample_distances_A_sqrt.dist", 
                       delim = "\t", col_names = FALSE)
#load metadata

# View first rows
head(metadata)

# Extract sample names (column 1)
sample_names <- dist_raw[[1]]

# Extract only the distance values: skip columns 1 (sample name) and 2 (ID)
dist_values <- as.matrix(dist_raw[, -(1:2)])

# Assign row and column names
rownames(dist_values) <- sample_names
colnames(dist_values) <- sample_names  # Assumes symmetric distances in same order

# Confirm it's square
print(dim(dist_values))  # Should be 32 x 32

# Create the distance object
dist_object <- as.dist(dist_values)


# Now convert to distance object
dist_object <- as.dist(dist_values)

metadata <- metadata %>% subset(Sample != 'HBL4')

# Run PERMANOVA
adonis_result <- adonis2(dist_object ~ Location, data = metadata)
print(adonis_result)

###Morphotype instead of island
# Run PERMANOVA by Morphotype
adonis_result_morphotype <- adonis2(dist_object ~ Morph, data = metadata)
print(adonis_result_morphotype)

adonis2(dist_object ~ Morph * Location, data = metadata)


# Load the necessary libraries
#library(pairwiseAdonis)
#install.packages('devtools')
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Perform pairwise PERMANOVA
pairwise_results <- pairwise.adonis2(dist_object ~ Location:Morph, data = metadata)
print(pairwise_results)


