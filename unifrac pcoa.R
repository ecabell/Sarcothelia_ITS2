# Load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ape)  # for pcoa
library(tibble)
library(ggrepel)
# ---------------------------
# Step 1: Read UniFrac distance matrix
# ---------------------------
dist_raw <- read_delim("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist", 
                       delim = "\t", col_names = FALSE)

sample_ids <- dist_raw$X1
dist_values <- dist_raw[, -(1:2)]

dist_matrix <- as.matrix(dist_values)
rownames(dist_matrix) <- sample_ids
colnames(dist_matrix) <- sample_ids
dist_matrix <- as.dist(dist_matrix)

# ---------------------------
# Step 2: Perform PCoA
# ---------------------------
pcoa_result <- pcoa(dist_matrix)

# ---------------------------
# Step 3: Prepare data for plotting
# ---------------------------
pcoa_df <- as_tibble(pcoa_result$vectors) %>%
  mutate(SampleID = rownames(pcoa_result$vectors))

# ---------------------------
# Step 4: Plot the first two axes
# ---------------------------
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, label = SampleID)) +
  geom_point(size = 3) +
  geom_text_repel(vjust = -1, size = 3) +
  labs(title = expression("PCoA of "*italic("Symbiodinium")*" (clade A) UniFrac Distance Matrix"),
       x = paste0("PCoA 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_minimal()

mean_dist <- mean(pcoa_df$Distance)
sd_dist <- sd(pcoa_df$Distance)
pcoa_df <- pcoa_df %>%
  mutate(Significance = ifelse(Distance > mean_dist + 2 * sd_dist, "Far", "Close"))

center <- colMeans(pcoa_df[, c("Axis.1", "Axis.2")])
cov_matrix <- cov(pcoa_df[, c("Axis.1", "Axis.2")])
pcoa_df$Mahalanobis <- mahalanobis(pcoa_df[, c("Axis.1", "Axis.2")], center, cov_matrix)

# Use a chi-square distribution threshold for 2 dimensions (e.g., 95% cutoff)
threshold <- qchisq(0.95, df = 2)
pcoa_df <- pcoa_df %>%
  mutate(Significance = ifelse(Mahalanobis > threshold, "Far", "Close"))



# Hierarchical clustering on the UniFrac distance matrix
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram
plot(hc, labels = sample_ids, main = "Hierarchical Clustering of Samples")

library(vegan)

# Assume you have metadata with a column 'Group' (e.g., treatment, site, etc.)
metadata <- read_csv("~/Documents/Zoox PCR/ITS2_Profile_Metadata.csv")  # must match sample_ids
metadata <- metadata %>% rename(SampleID = ITS2_Profile)
metadata <- metadata %>% filter(SampleID %in% sample_ids)

# Run PERMANOVA
permanova <- adonis2(dist_matrix ~ Morph, data = metadata)
print(permanova)
# Run PERMANOVA
permanova <- adonis2(dist_matrix ~ Location, data = metadata)
print(permanova)

# Test homogeneity of group variances (dispersion)
disp <- betadisper(dist_matrix, metadata$Morph)

# Plot and test
plot(disp)
boxplot(disp, main = "Dispersion by Morph")
anova(disp)  # ANOVA test for differences in dispersion

disp <- betadisper(dist_matrix, metadata$Location)

# Plot and test
plot(disp)
boxplot(disp, main = "Dispersion by Island")
anova(disp)  # ANOVA test for differences in dispersion

pcoa_df$Cluster <- factor(clusters)
pcoa_df$Group <- metadata$Group[match(pcoa_df$SampleID, metadata$SampleID)]

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  labs(title = expression("PCoA of "*italic("UniFrac")*" Distance Matrix")) +
  theme_minimal()





######### between samples

# Load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ape)  # for pcoa
library(tibble)
library(ggrepel)
# ---------------------------
# Step 1: Read UniFrac distance matrix
# ---------------------------
dist_raw <- read_delim("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_sample_distances/A/20250605T195712_unifrac_sample_distances_A_sqrt.dist", 
                       delim = "\t", col_names = FALSE)

sample_ids <- dist_raw$X1
dist_values <- dist_raw[, -(1:2)]

dist_matrix <- as.matrix(dist_values)
rownames(dist_matrix) <- sample_ids
colnames(dist_matrix) <- sample_ids
dist_matrix <- as.dist(dist_matrix)

# ---------------------------
# Step 2: Perform PCoA
# ---------------------------
pcoa_result <- pcoa(dist_matrix)

# ---------------------------
# Step 3: Prepare data for plotting
# ---------------------------
pcoa_df <- as_tibble(pcoa_result$vectors) %>%
  mutate(SampleID = rownames(pcoa_result$vectors))

# ---------------------------
# Step 4: Plot the first two axes
# ---------------------------
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, label = SampleID)) +
  geom_point(size = 3) +
  geom_text_repel(vjust = -1, size = 3) +
  labs(title = expression("PCoA of Samples UniFrac Distance Matrix"),
       x = paste0("PCoA 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_minimal()

center <- colMeans(pcoa_df[, c("Axis.1", "Axis.2")])
cov_matrix <- cov(pcoa_df[, c("Axis.1", "Axis.2")])
pcoa_df$Mahalanobis <- mahalanobis(pcoa_df[, c("Axis.1", "Axis.2")], center, cov_matrix)

# Use a chi-square distribution threshold for 2 dimensions (e.g., 95% cutoff)
threshold <- qchisq(0.95, df = 2)
pcoa_df <- pcoa_df %>%
  mutate(Significance = ifelse(Mahalanobis > threshold, "Far", "Close"))

# Hierarchical clustering on the UniFrac distance matrix
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram
plot(hc, labels = sample_ids, main = "Hierarchical Clustering of Samples")
clusters <- cutree(hc, k = 3)

# Test homogeneity of group variances (dispersion)
#load metadata
metadata <- data.frame(SampleID = sample_names) %>%
  mutate(
    Location = case_when(
      grepl("K", SampleID) ~ "Kauai",
      grepl("O", SampleID) ~ "Oahu",
      grepl("H", SampleID) ~ "Big Island",
      TRUE ~ NA_character_
    ),
    Morph = case_when(
      grepl("BR", SampleID) ~ "Brown",
      grepl("BL", SampleID) ~ "Blue",
      TRUE ~ NA_character_
    )
  )

disp <- betadisper(dist_matrix, metadata$Morph)

# Plot and test
plot(disp)
boxplot(disp, main = "Dispersion by Morph")
anova(disp)  # ANOVA test for differences in dispersion

disp <- betadisper(dist_matrix, metadata$Location)

# Plot and test
plot(disp)
boxplot(disp, main = "Dispersion by Island")
anova(disp)  # ANOVA test for differences in dispersion

pcoa_df$Cluster <- factor(clusters)
pcoa_df$Group <- metadata$Group[match(pcoa_df$SampleID, metadata$SampleID)]

ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Significance, label = SampleID)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('black', 'red')) + 
  geom_text(vjust = -1, size = 3) +
  labs(title = expression("PCoA of "*italic("UniFrac")*" Distance Matrix")) +
  theme_minimal()

