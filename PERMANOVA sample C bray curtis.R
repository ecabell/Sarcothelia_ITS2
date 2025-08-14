# Load required packages
library(vegan)

# Load distance matrix
dist_data <- read.table("~/Documents/SEdmondsoni_Zoox_Seq/20250422T005935/between_sample_distances/C/20250422T005935_braycurtis_sample_distances_C_sqrt.dist", header = FALSE, row.names = 1, fill = TRUE)

# Convert to a full distance matrix
sample_names <- rownames(dist_data)
dist_matrix <- as.dist(as.matrix(dist_data[, -1]))

# Create metadata dataframe
metadata <- data.frame(Sample = sample_names)
metadata$Island <- ifelse(grepl("K", metadata$Sample), "Kauai",
                          ifelse(grepl("O", metadata$Sample), "Oahu",
                                 ifelse(grepl("H", metadata$Sample), "Big Island", NA)))

# Run PERMANOVA
adonis_result <- adonis2(dist_matrix ~ Island, data = metadata)
print(adonis_result)

# Load distance matrix #no sqrt
dist_data <- read.table("~/Documents/SEdmondsoni_Zoox_Seq/20250422T005935/between_sample_distances/C/20250422T005935_braycurtis_sample_distances_C_no_sqrt.dist", header = FALSE, row.names = 1, fill = TRUE)

# Convert to a full distance matrix
sample_names <- rownames(dist_data)
dist_matrix <- as.dist(as.matrix(dist_data[, -1]))

# Create metadata dataframe
metadata <- data.frame(Sample = sample_names)
metadata$Island <- ifelse(grepl("K", metadata$Sample), "Kauai",
                          ifelse(grepl("O", metadata$Sample), "Oahu",
                                 ifelse(grepl("H", metadata$Sample), "Big Island", NA)))

# Run PERMANOVA
adonis_result <- adonis2(dist_matrix ~ Island, data = metadata)
print(adonis_result)

###Morphotype instead of island
# Create metadata
metadata <- data.frame(Sample = sample_names)
metadata$Morphotype <- ifelse(grepl("BL", metadata$Sample), "Blue", "Brown")

# Run PERMANOVA by Morphotype
adonis_result_morphotype <- adonis2(dist_matrix ~ Morphotype, data = metadata)
print(adonis_result_morphotype)

