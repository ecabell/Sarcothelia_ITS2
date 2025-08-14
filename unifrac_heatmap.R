# Styled UniFrac Distance Heatmap (matches ML figure style)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Load the same UniFrac distance matrix ---
data <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist", 
                   sep = "\t", header = FALSE, row.names = 1)
dist_matrix <- as.matrix(data[, -1])
sample_names <- rownames(dist_matrix)
n_samples <- nrow(dist_matrix)

# Symmetrize the matrix
symmetric_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
rownames(symmetric_matrix) <- colnames(symmetric_matrix) <- sample_names
for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    symmetric_matrix[i, j] <- dist_matrix[i, j]
    symmetric_matrix[j, i] <- dist_matrix[i, j]
  }
}

# Round display values
display_matrix <- round(symmetric_matrix, 3)

# Optional: apply sample label mapping
# label_map <- c("1234" = "A3", "5678" = "A4")
# rownames(display_matrix) <- colnames(display_matrix) <- sapply(rownames(display_matrix), function(x) {
#   if (!is.null(label_map[[x]])) label_map[[x]] else x
# })

# --- Define color palette and annotation style ---
col_fun <- colorRamp2(c(0, max(display_matrix)), c("seashell", "#08306b"))
ht_opt(HEATMAP_COLUMN_TITLE_PADDING = unit(2, "mm"), ADD = TRUE)

# Build heatmap object
# Build heatmap object with hierarchical clustering
ht <- Heatmap(symmetric_matrix,
              name = "UniFrac",
              col = col_fun,
              cluster_rows = TRUE,         # ← enable clustering
              cluster_columns = TRUE,      # ← enable clustering
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 11, fontfamily = "Times"),
              row_names_gp = gpar(fontsize = 11, fontfamily = "Times"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- display_matrix[i, j]
                font_col <- if (val >= 0.04) "white" else "black"
                font_weight <- if (val >= 0.04) "bold" else "plain"
                grid.text(sprintf("%.3f", val), x, y,
                          gp = gpar(fontsize = 10, col = font_col, fontface = font_weight, fontfamily = "Times"))
              },
              column_title = "B) UniFrac Distance Matrix Heatmap (Clustered)",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times"),
              width = unit(13, "cm"),
              height = unit(13, "cm"))

# --- Save as PDF and PNG with Times font ---
pdf("unifrac_heatmap_times.pdf", family = "Times", width = 8, height = 8)
draw(ht, heatmap_legend_side = "right", padding = unit(c(10, 5, 10, 10), "mm"))
dev.off()

png("unifrac_heatmap_times.png", width = 2400, height = 2400, res = 300, family = "Times")
draw(ht, heatmap_legend_side = "right")
dev.off()

cat("Heatmap saved as:\n")
cat(" - unifrac_heatmap_times.pdf\n")
cat(" - unifrac_heatmap_times.png\n")
