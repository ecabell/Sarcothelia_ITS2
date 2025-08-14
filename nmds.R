
# Load required libraries
library(readr)
library(scico)
library(vegan)

# --- 1. Load Distance Matrix ---
dist_raw <- read.delim("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_sample_distances/A/20250605T195712_braycurtis_sample_distances_A_sqrt.dist", 
                       header = FALSE)
sample_names <- dist_raw[, 1]
dist_values <- dist_raw[, -(1:2)]
dist_matrix <- as.matrix(dist_values)
rownames(dist_matrix) <- sample_names
colnames(dist_matrix) <- sample_names
dist_object <- as.dist(dist_matrix)

# --- 2. Create Metadata Data Frame ---
metadata <- data.frame(SampleID = sample_names, stringsAsFactors = FALSE) %>%
  mutate(
    Island = case_when(
      grepl("^K", SampleID) ~ "Kauai",
      grepl("^O", SampleID) ~ "Oahu",
      grepl("^H", SampleID) ~ "Big Island",
      TRUE ~ NA_character_
    ),
    Morphotype = case_when(
      grepl("BL", SampleID) ~ "Blue",
      grepl("BR", SampleID) ~ "Brown",
      TRUE ~ NA_character_
    )
  )
# Ensure factors
metadata$Island <- factor(metadata$Island)
metadata$Morphotype <- factor(metadata$Morphotype)

# --- 3. Run NMDS ---
set.seed(123)
nmds_result <- metaMDS(dist_object, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
cat("Stress =", nmds_result$stress, "\n")

# --- 4. Extract Site Scores & Merge Metadata ---
site_scores <- as.data.frame(scores(nmds_result, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
site_scores <- merge(site_scores, metadata, by = "SampleID")

# --- 5. Plot NMDS: Color by Island, Shape by Morphotype ---
island_colors <- setNames(scico(length(levels(metadata$Island)), palette = "hawaii"), levels(metadata$Island))
morph_shapes <- setNames(c(19, 17), levels(metadata$Morphotype))

plot(nmds_result, type = "n",
     main = "NMDS: Island (color) & Morphotype (shape)")
points(site_scores$NMDS1, site_scores$NMDS2,
       col = island_colors[site_scores$Island],
       pch = morph_shapes[site_scores$Morphotype])
# Draw ellipses by island
for (isl in levels(site_scores$Island)) {
  ordiellipse(nmds_result,
              groups = site_scores$Island,
              kind = "sd", conf = 0.95,
              show.groups = isl,
              col = island_colors[isl], lwd = 2)
}
legend("topright", legend = levels(metadata$Island),
       col = island_colors, pch = 19, title = "Island")
legend("bottomright", legend = levels(metadata$Morphotype),
       pch = morph_shapes, title = "Morphotype")

# --- 6. PERMANOVA Analysis ---
adonis_result <- adonis2(dist_object ~ Island + Morphotype, data = metadata, permutations = 999)
print(adonis_result)
# --- 5. Plot NMDS by Group (e.g., Morphotype) ---
group_factor <- as.factor(metadata$Morphotype)
colors <- rainbow(length(unique(group_factor)))[group_factor]

plot(nmds_result, type = "n", main = "NMDS by Morphotype")
points(nmds_result, display = "sites", col = colors, pch = 19)
ordiellipse(nmds_result, group_factor, kind = "sd", conf = 0.95)

legend("bottomleft", legend = levels(group_factor), col = rainbow(length(unique(group_factor))), pch = 19)
par(xpd = TRUE)
ordiellipse(nmds_result, group_factor, kind = "sd", conf = 0.95)
par(xpd = FALSE)

# Get site scores
site_scores <- scores(nmds_result, display = "sites")

# Determine limits
xlims <- range(site_scores[,1]) + c(-0.1, 0.1)
ylims <- range(site_scores[,2]) + c(-0.1, 0.1)

# Replot with custom limits
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS by Morphotype")
points(nmds_result, display = "sites", col = colors, pch = 19)
ordiellipse(nmds_result, group_factor, kind = "sd", conf = 0.95)
legend("bottomleft", legend = levels(group_factor), col = rainbow(length(unique(group_factor))), pch = 19)

# Define custom colors
morph_colors <- c("Blue" = "blue", "Brown" = "brown")

# Extract Morphotype as a factor
group_factor <- as.factor(metadata$Morphotype)

# Match color vector to morphotype
point_colors <- morph_colors[as.character(group_factor)]

# Plot with manual limits if needed
site_scores <- scores(nmds_result, display = "sites")
xlims <- range(site_scores[,1]) + c(-0.1, 0.1)
ylims <- range(site_scores[,2]) + c(-0.1, 0.1)

# Plot base
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS by Morphotype")
points(nmds_result, display = "sites", col = point_colors, pch = 19)

# Add ellipses manually per group
for (morph in levels(group_factor)) {
  ordiellipse(nmds_result, group_factor,
              kind = "sd", conf = 0.95, show.groups = morph,
              col = morph_colors[morph], lwd = 2)
}

# Add legend
legend("bottomleft", legend = names(morph_colors), col = morph_colors, pch = 19)


# Plot with manual limits if needed
site_scores <- scores(nmds_result, display = "sites")
xlims <- range(site_scores[,1]) + c(-0.1, 0.3)
ylims <- range(site_scores[,2]) + c(-0.1, 0.1)
# Factor and color map
group_factor <- as.factor(metadata$Island)
island_levels <- levels(group_factor)
island_colors <- setNames(rainbow(length(island_levels)), island_levels)

# Color points
point_colors <- island_colors[as.character(group_factor)]

# Plot
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS by Island")
points(nmds_result, display = "sites", col = point_colors, pch = 19)

# Ellipses colored by island
for (grp in island_levels) {
  ordiellipse(nmds_result, group_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = island_colors[grp], lwd = 2)
}

# Legend
legend("bottomleft", legend = island_levels, col = island_colors, pch = 19)

# Factor levels
island_factor <- as.factor(metadata$Island)
morph_factor <- as.factor(metadata$Morphotype)

library(viridisLite)
library(scico)
# Define color map for islands
island_levels <- levels(island_factor)
#island_colors <- setNames(scico(length(island_levels), palette  = "berlin"), island_levels)
#hawaii #vanimo #also goood palettes #batlow

island_colors <- setNames(viridis(length(island_levels), option = "inferno"), island_levels)

# Replace third color with a stronger one (e.g., "red" or "darkorange")
island_colors[3] <- "darkgreen"
point_colors <- island_colors[as.character(island_factor)]
island_shapes <- setNames(c(19, 17, 15), island_levels)

# Define shape map for morphotype
morph_levels <- levels(morph_factor)
morph_shapes <- setNames(c(19, 17), morph_levels)  # e.g., 19 = solid circle, 17 = triangle
point_shapes <- morph_shapes[as.character(morph_factor)]
# Define shape map for islands
point_shapes2 <- island_shapes[as.character(island_factor)]
point_colors2 <- morph_colors[as.character(morph_factor)]

# Plot limits
site_scores <- scores(nmds_result, display = "sites")
xlims <- range(site_scores[,1]) + c(-0.5, 0.3)
ylims <- range(site_scores[,2]) + c(-0.1, 0.1)

# Plot
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS: Color by Island, Shape by Morphotype")
points(nmds_result, display = "sites", col = point_colors, pch = point_shapes)

# Add ellipses by island
for (grp in island_levels) {
  ordiellipse(nmds_result, island_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = island_colors[grp], lwd = 2)
}

# Add legends
legend("bottomleft", inset = 0.05, legend = island_levels, col = island_colors, pch = 19, title = "Island")
legend("bottomright", inset = 0.05, legend = morph_levels, pch = unlist(morph_shapes), title = "Morphotype")

# Plot
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS: Color by Morphotype, Shape by Island")
points(nmds_result, display = "sites", col = point_colors2, pch = point_shapes2)

# Add ellipses by morph
for (grp in morph_levels) {
  ordiellipse(nmds_result, morph_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = morph_colors[grp], lwd = 2)
}

# Add ellipses by island
for (grp in island_levels) {
  ordiellipse(nmds_result, island_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = island_colors[grp], lwd = 1, lty = "dotted")
}

# Assume you've already computed centroids_island and have island_colors
points(centroids_island$NMDS1, centroids_island$NMDS2, pch = 4, col = island_colors[centroids_island$Island], cex = 1.5)

text(centroids_island$NMDS1,
     centroids_island$NMDS2,
     labels = centroids_island$Island,
     col = island_colors[centroids_island$Island],  # color match
     cex = 0.75,                                     # smaller text
     pos = 2)                                        # next to the point

# Add legends
legend("topleft", inset = 0.05, legend = morph_levels, col = morph_colors, pch = 19, title = "Morphotype")
legend("topright", inset = 0.05, legend = island_levels, col = island_colors, pch = unlist(island_shapes), title = "Island")

# Plot
plot(nmds_result, type = "n", xlim = xlims, ylim = ylims, main = "NMDS: Color by Island, Shape by Morphotype")
points(nmds_result, display = "sites", col = point_colors, pch = point_shapes)

# Add ellipses by island
for (grp in island_levels) {
  ordiellipse(nmds_result, island_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = island_colors[grp], lwd = 2)
}

# Add ellipses by morph
for (grp in morph_levels) {
  ordiellipse(nmds_result, morph_factor, kind = "sd", conf = 0.95,
              show.groups = grp, col = morph_colors[grp], lwd = 1, lty = 'dotted')
}

# Add legends
legend("topright", inset = 0.05, legend = island_levels, col = island_colors, pch = 19, title = "Island")
legend("bottomright", inset = 0.05, legend = morph_levels, pch = unlist(morph_shapes), title = "Morphotype")


# Ensure Island and Morphotype are factors
metadata$Island <- as.factor(metadata$Island)
metadata$Morphotype <- as.factor(metadata$Morphotype)

# Run PERMANOVA
adonis_result <- adonis2(dist_object ~ Island + Morphotype, data = metadata, permutations = 999)

# View results
print(adonis_result)

adonis2(dist_object ~ Island, data = metadata)
adonis2(dist_object ~ Morphotype, data = metadata)

