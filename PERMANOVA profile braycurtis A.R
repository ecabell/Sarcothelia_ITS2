# Load libraries
library(readr)
library(dplyr)
library(vegan)

# -----------------------
# Step 1: Read distance matrix
# -----------------------

# Read tab-delimited .dist file (no header, first column = profile, second = index)
dist_raw <- read_delim("~/Documents/Docs/Zoox PCR/SEdmondsoni_Zoox_Seq/20250422T005935/between_profile_distances/A/20250422T005935_braycurtis_profile_distances_A_sqrt.dist", delim = "\t", col_names = FALSE)

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
metadata <- read_csv("~/Documents/Docs/Zoox PCR/SEdmondsoni_Zoox_Seq/20250422T005935/Grouped_ITS2_Profile_Data.csv")

# Collapse to one row per ITS2 profile using highest abundance
metadata_unique <- metadata %>%
  group_by(ITS2_Profile) %>%
  slice_max(Value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Set rownames (base R, avoids tibble warnings later)
metadata_unique_df <- as.data.frame(metadata_unique)
rownames(metadata_unique_df) <- metadata_unique_df$ITS2_Profile
metadata_unique_df$ITS2_Profile <- NULL

# -----------------------
# Step 3: Align metadata and distance matrix
# -----------------------

# Get common profile names
common_names <- intersect(rownames(metadata_unique_df), labels(dist_matrix))

# Subset metadata and distance matrix
metadata_final <- metadata_unique_df[common_names, ]
dist_matrix_final <- as.dist(as.matrix(dist_matrix)[common_names, common_names])

# -----------------------
# Step 4: Run PERMANOVA
# -----------------------

# Test effect of Island
adonis_island <- adonis2(dist_matrix_final ~ Island, data = metadata_final)
print(adonis_island)

# Test effect of Morphotype
adonis_morph <- adonis2(dist_matrix_final ~ Morphotype, data = metadata_final)
print(adonis_morph)

