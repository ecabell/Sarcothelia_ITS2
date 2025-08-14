# Load required packages
library(vegan)

# Load distance matrix
dist_data <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_sample_distances/A/20250605T195712_unifrac_sample_distances_A_sqrt.dist", header = FALSE, row.names = 1, fill = TRUE)

# Convert to a full distance matrix
sample_names <- rownames(dist_data)
dist_matrix <- as.dist(as.matrix(dist_data[, -1]))

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

# Run PERMANOVA
adonis_result <- adonis2(dist_matrix ~ Location, data = metadata)
print(adonis_result)

## Summary:
#While your PCoA shows visual clustering, the differences in microbial communities between islands are not statistically significant at the p < 0.05 level.

###Morphotype instead of island
# Run PERMANOVA by Morphotype
adonis_result_morphotype <- adonis2(dist_matrix ~ Morph, data = metadata)
print(adonis_result_morphotype)
 

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
