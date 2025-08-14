# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# ---------------------------
# Step 1: Read UniFrac distance matrix
# ---------------------------
dist_raw <- read_delim("~/Documents/SEdmondsoni_Zoox_Seq/20250422T005935/between_sample_distances/C/20250422T005935_braycurtis_sample_distances_C_sqrt.dist", 
                       delim = "\t", col_names = FALSE)

sample_ids <- dist_raw$X1
dist_values <- dist_raw[, -(1:2)]

dist_matrix <- as.matrix(dist_values)
rownames(dist_matrix) <- sample_ids
colnames(dist_matrix) <- sample_ids
dist_matrix <- as.dist(dist_matrix)

# ---------------------------
# Step 2: Run PCoA
# ---------------------------
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
pcoa_coords$Sample_Name <- rownames(pcoa_coords)

# ---------------------------
# Step 3: Load and clean metadata
# ---------------------------
# Read the sample-level metadata file (UTF-16 and tab-separated)
meta_raw <- read_delim("~/Documents/SEdmondsoni_Zoox_Seq/20250422T005935/TypebySample.csv", 
                       delim = "\t", col_names = FALSE, locale = locale(encoding = "UTF-16"))

# Split and rename columns
colnames(meta_raw) <- c("Sample_Name", "ITS2_Profile", "Value")

# Add Island and Morphotype from sample names
meta_clean <- meta_raw %>%
  mutate(
    Island = case_when(
      grepl("^K", Sample_Name) ~ "Kauai",
      grepl("^H", Sample_Name) ~ "Big Island",
      grepl("^O", Sample_Name) ~ "Oahu",
      TRUE ~ "Unknown"
    ),
    Morphotype = case_when(
      grepl("BL", Sample_Name) ~ "Blue",
      grepl("BR", Sample_Name) ~ "Brown",
      TRUE ~ "Unknown"
    )
  )

# ---------------------------
# Step 4: Merge with PCoA coordinates
# ---------------------------
plot_data <- left_join(pcoa_coords, meta_clean, by = "Sample_Name")

# ---------------------------
# Step 5: Plot PCoA
# ---------------------------
var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)

# Calculate centroids for ellipse labels
centroids_morph <- plot_data %>%
  group_by(Morphotype) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

centroids_island <- plot_data %>%
  group_by(Island) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop") %>%
  mutate(label = paste(Island, "centroid"))


# Plot with labeled ellipses
ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  # Morphotype ellipses (dashed)
  stat_ellipse(aes(group = Morphotype, color = Morphotype), linetype = "dashed", size = 0.7) +
  # Island ellipses (solid black)
 # stat_ellipse(aes(group = Island), color = "black", linetype = "solid", size = 0.7, alpha = 0.4) +
  # Points
  geom_point(aes(color = Morphotype, shape = Island), size = 2.5, stroke = 1.2) +
  # Ellipse labels
  geom_text(data = centroids_morph, aes(label = Morphotype, color = Morphotype), fontface = "bold", size = 3, hjust = -.5) +
  geom_text(data = centroids_island, aes(label = Island), color = "black", fontface = "italic", size = 3, vjust = 1.5, hjust = 1) +
  scale_color_manual(values = c("Blue" = "turquoise3", "Brown" = "coral3")) +
  theme_minimal() +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    title = "PCoA of Samples (Bray-Curtis)\nColored by Morphotype, Shaped by Island",
    color = "Morphotype",
    shape = "Island"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )

#That warning from stat_ellipse() means that at least one group 
#has too few (or no) points, or it contains NA values — 
#which causes the internal ellipse-fitting function to fail.
