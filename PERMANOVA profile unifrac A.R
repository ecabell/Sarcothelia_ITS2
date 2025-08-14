-# Load libraries
library(readr)
library(dplyr)
library(vegan)

# -----------------------
# Step 1: Read distance matrix
# -----------------------

# Read tab-delimited .dist file (no header, first column = profile, second = index)
dist_raw <- read_delim("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist", delim = "\t", col_names = FALSE)

# Extract profile names
profile_names <- dist_raw$X1

# Remove index column and profile name column
dist_values <- dist_raw[, -(1:2)]

# Convert to matrix and assign row/column names
dist_matrix <- as.matrix(dist_values)
rownames(dist_matrix) <- profile_names
colnames(dist_matrix) <- profile_names

# Convert to dist object
dist_matrix <- as.dist(dist_matrix)

# -----------------------
# Step 2: Read and prepare metadata
# -----------------------

# Read metadata
# Create metadata dataframe
metadata <- data.frame(ProfileName = profile_names) %>%
  mutate(
    Type = ifelse(grepl("A4", ProfileName), "A4", "A3")
  )
# -----------------------
# Step 4: Run PERMANOVA
# -----------------------

# Test effect of Type
adonis_island <- adonis2(dist_matrix ~ Type, data = metadata)
print(adonis_island)

