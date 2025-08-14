# Load required packages
library(ape)
library(pheatmap)

# Set seed for reproducibility
set.seed(123)

# --- Step 1: Read aligned FASTA file ---
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# --- Step 2: Compute uncorrected p-distance matrix ---
dist_p <- dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)

# --- Step 3: Build unrooted NJ tree ---
nj_tree <- unroot(nj(dist_p))

# --- Step 4: Bootstrap the NJ tree (100 replicates) ---
boot_nj <- boot.phylo(nj_tree, alignment, FUN = function(x) {
  tree <- nj(dist.dna(x, model = "raw"))
  unroot(multi2di(tree))
}, B = 100)

# --- Step 5: Rename tip labels ---
label_map <- c(
  "699985" = "A4ac", "760355" = "A4eo", "33881" = "A3bt", "22415" = "A3bb",
  "22385" = "A3", "34175" = "A3fa", "34422" = "A3u", "36733" = "A3ah",
  "35231" = "A3ad", "502153" = "A3ez", "363640" = "A3ey", "363655" = "A3fb",
  "22469" = "A3g"
)
nj_tree$tip.label <- sapply(nj_tree$tip.label, function(x) {
  if (!is.null(label_map[[x]])) label_map[[x]] else x
})

# --- Step 6: Plot the NJ tree with bootstrap values ---
ntips <- Ntip(nj_tree)
par(mar = c(6, 4, 4, 2))
plot(nj_tree, show.tip.label = TRUE, cex = 0.8, main = "A) Neighbor-Joining Tree with Bootstrap Support ≥50%", y.lim = c(-2, ntips))
nodelabels(text = ifelse(boot_nj >= 50, boot_nj, ""), frame = "rect", cex = 0.7, font = 2)

# --- Step 7: Add scale bar for uncorrected p-distance ---
usr <- par("usr")
x_pos <- usr[1] + 0.003
y_pos <- usr[3] + 1.2
add.scale.bar(x = x_pos, y = y_pos + 1, length = 0.06, lwd = 1, cex = 0.8)
text(x = x_pos + 0.03, y = y_pos + 0.7,
     labels = "0.01 substitutions/site (uncorrected p-distance)", cex = 0.9)

###### Heatmap
dist_matrix <- as.matrix(dist_p)

rownames(dist_matrix) <- colnames(dist_matrix) <- sapply(rownames(dist_matrix), function(x) {
  if (!is.null(label_map[[x]])) label_map[[x]] else x
})

# --- Step 8: Plot the heatmap ---
pheatmap(dist_matrix,
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         fontsize_number = 8,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "B) ITS2 Type Uncorrected p-distance Heatmap", 
         legend_position = "right")

library(ComplexHeatmap)
library(circlize)

# Prepare the matrix with updated labels and rounded values
dist_matrix <- as.matrix(dist_p)
rownames(dist_matrix) <- colnames(dist_matrix) <- sapply(rownames(dist_matrix), function(x) {
  if (!is.null(label_map[[x]])) label_map[[x]] else x
})
display_matrix <- round(dist_matrix, 3)

# Define color palette
col_fun <- colorRamp2(c(0, max(dist_matrix)), c("seashell", "#08306b"))
# Set global spacing options to reduce white box size
ht_opt(HEATMAP_COLUMN_TITLE_PADDING = unit(2, "mm"), ADD = TRUE)

cell_fun = function(j, i, x, y, width, height, fill) {
  value <- display_matrix[i, j]
  font_col <- if (value >= 0.04) "ivory" else "black"
  font_weight <- if (value >= 0.04) "bold" else "plain"
  
  grid.text(sprintf("%.3f", value), x, y,
            gp = gpar(fontsize = 10,
                      col = font_col,
                      fontface = font_weight,
                      fontfamily = "Times"))
}
# Now create the heatmap
# Increase cell size and slightly reduce font size for better legibility
ht <- Heatmap(dist_matrix,
              name = "p-distance",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_names = TRUE,
              show_row_names = TRUE,
              column_names_gp = gpar(fontsize = 11, fontfamily = 'Times'),
              row_names_gp = gpar(fontsize = 11, fontfamily = 'Times'),
              cell_fun = function(j, i, x, y, width, height, fill) {
                value <- display_matrix[i, j]
                font_col <- if (value >= 0.04) "white" else "black"
                grid.text(sprintf("%.3f", value), x, y,                           
                          gp = gpar(fontsize = 10, col = font_col, fontfamily = "Times"))

              },
              column_title = "B) ITS2 Type Uncorrected p-distance Heatmap",
              column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = 'Times'),
              width = unit(16, "cm"),   # Slightly wider
              height = unit(14, "cm"))  # Slightly taller

draw(ht)

pdf("heatmap_times.pdf", family = "Times", width = 8, height = 8)
draw(ht)
dev.off()

png("heatmap_times.png", width = 2400, height = 2400, res = 300, family = "Times")
draw(ht)
dev.off()
