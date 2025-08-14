# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# ---------------------------
# Step 1: Read Bray Curtis distance matrix
# ---------------------------
dist_raw <- read_delim("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_sample_distances/A/20250605T195712_braycurtis_sample_distances_A_sqrt.dist", 
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
#load metadata (metadata r file)

# ---------------------------
# Step 4: Merge with PCoA coordinates
# ---------------------------


plot_data <- left_join(pcoa_coords, metadata, by = 'Sample_Name')

# ---------------------------
# Step 5: Plot PCoA
# ---------------------------
var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)

# Calculate centroids for ellipse labels
centroids_Morph <- plot_data %>%
  group_by(Morph) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

centroids_Location <- plot_data %>%
  group_by(Location) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop") %>%
  mutate(label = paste(Location, "centroid"))


# Plot with labeled ellipses
ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  # Morph ellipses (dashed)
  stat_ellipse(aes(group = Morph, color = Morph), linetype = "dashed", size = 0.7) +
  # Location ellipses (solid black)
  stat_ellipse(aes(group = Location), color = "black", linetype = "solid", size = 0.7, alpha = 0.4) +
  # Points
  geom_point(aes(color = , shape = Location), size = 2.5, stroke = 1.2) +
  # Ellipse labels
  geom_text(data = centroids_Morph, aes(label = Morph, color = Morph), fontface = "bold", size = 3, vjust = -1, hjust = -.5, show.legend = FALSE) +
  geom_text(data = centroids_Location, aes(label = Location), color = "black", fontface = "italic", size = 3, vjust = 1.5, show.legend =  FALSE) +
  scale_color_manual(values = c("Blue" = "blue", "Brown" = "brown")) +
  theme_minimal() +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    title = "PCoA of Samples (Bray-Curtis)\nColored by Morph, Shaped by Location",
    color = "Morph",
    shape = "Location"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )

theme_nmds_style <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  # Morph ellipses (dashed)
  stat_ellipse(aes(group = Morph, color = Morph), linetype = "dashed", size = 0.5) +
  # Location ellipses (solid black)
  stat_ellipse(aes(group = Location), color = "black", linetype = "solid", size = 0.5, alpha = 0.4) +
  # Points
  geom_point(aes(color = Morph, shape = Location), size = 2.5, stroke = 1.2) +
  # Ellipse labels
  geom_text(data = centroids_Morph, aes(label = Morph, color = Morph), fontface = "bold", size = 3, vjust = -1, hjust = -.5, show.legend = FALSE) +
  geom_text(data = centroids_Location, aes(label = Location), color = "black", fontface = "italic", size = 3, vjust = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("Blue" = "blue", "Brown" = "brown")) +
  theme_nmds_style +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    title = "PCoA: Color by Morph, Shape by Location", 
    color = "Morph",
    shape = "Location") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 12),
    legend.position = "right")


#####color by Location shape by Morpho
# Calculate centroids for updated legend positions
centroids_Location <- plot_data %>%
  group_by(Location) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

centroids_Morph <- plot_data %>%
  group_by(Morph) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")

# Define custom colors for Locations (no blue)
Location_colors <- c(
  "Kauai" = "darkorange2",
  "Big Location" = "forestgreen",
  "Oahu" = "purple3"
)

ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  # Location ellipses (colored)
  stat_ellipse(aes(group = Location, color = Location), linetype = "solid", size = 0.7) +
  # Morph ellipses (dashed black)
  stat_ellipse(aes(group = Morph), color = "black", linetype = "dashed", size = 0.6, alpha = 0.4) +
  # Points: color by Location, shape by Morph
  geom_point(aes(color = Location, shape = Morph), size = 2.5, stroke = 1.1) +
  # Location centroid labels
  geom_text(data = centroids_Location, aes(label = Location, color = Location), fontface = "italic", size = 3, vjust = -1.5) +
  # Morph centroid labels
  geom_text(data = centroids_Morph, aes(label = Morph), color = "black", fontface = "bold", size = 3, vjust = 1.8) +
  scale_color_manual(values = Location_colors) +
  theme_minimal() +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    title = "PCoA of Samples (Bray-Curtis)\nColored by Location (No Blue), Shaped by Morph",
    color = "Location",
    shape = "Morph"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  )


pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
pcoa_result$eig


