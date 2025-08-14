# Styled Neighbor-Joining Tree with UniFrac Distance
# Author: Rewritten to match ML tree aesthetics

library(ape)

# === Load and prepare UniFrac distance matrix ===
unifrac_file <- "~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/20250605T195712_unifrac_profile_distances_A_sqrt.dist"
data <- read.table(unifrac_file, sep = "\t", header = FALSE, row.names = 1)

# Remove second column and extract distance matrix
dist_matrix <- as.matrix(data[, -1])
sample_names <- rownames(dist_matrix)
n_samples <- nrow(dist_matrix)

# Create symmetric matrix
symmetric_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
rownames(symmetric_matrix) <- sample_names
colnames(symmetric_matrix) <- sample_names
for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    symmetric_matrix[i, j] <- dist_matrix[i, j]
    symmetric_matrix[j, i] <- dist_matrix[i, j]
  }
}

# Convert to dist object and construct NJ tree
dist_obj <- as.dist(symmetric_matrix)
nj_tree <- nj(dist_obj)

# === Optional: Rename tip labels here if needed ===
# Example:
# label_map <- c("sample1" = "A1", "sample2" = "A2")
# nj_tree$tip.label <- sapply(nj_tree$tip.label, function(x) ifelse(!is.null(label_map[[x]]), label_map[[x]], x))

# === Plotting the styled NJ tree ===
pdf("unifrac_nj_tree_styled.pdf", width = 8, height = 7, family = "Times")
par(mar = c(6, 4, 6, 2), xpd = NA, family = "Times")

ntips <- Ntip(nj_tree)

# Main NJ tree plot
plot(unroot(nj_tree),
     show.tip.label = TRUE,
     cex = 1.2,
     font = 1,                   # Plain font
     main = "A) Neighbor-Joining Tree of ITS2 Types\n(Weighted UniFrac Distance)",
     y.lim = c(-2, ntips + 1))

# === Add UniFrac scale bar ===
usr <- par("usr")
scale_length <- 0.05  # Adjust to reflect meaningful UniFrac distance
x_start <- usr[1] + 0.02 * diff(usr[1:2])
y_start <- usr[3] + 2.7

segments(x_start, y_start, x_start + scale_length, y_start, lwd = 2)
text(x = x_start + scale_length / 2, y = y_start - 0.3,
     labels = "0.05 UniFrac distance", cex = 1.0, family = "Times")

dev.off()

# === Save tree and summary ===
write.tree(nj_tree, file = "unifrac_nj_tree.newick")

cat("\nTree Summary:\n")
print(nj_tree)
cat("\nNumber of tips:", Ntip(nj_tree), "\n")
cat("Tree is rooted:", is.rooted(nj_tree), "\n")
cat("Tree length:", sum(nj_tree$edge.length), "\n")

# === Optional: calculate R² of tree fit ===
cophenetic_dist <- cophenetic(nj_tree)
orig <- as.vector(symmetric_matrix)
coph <- as.vector(cophenetic_dist)
cor_val <- cor(orig, coph)
r2 <- cor_val^2
cat(sprintf("Variance explained (R²) by NJ tree: %.3f\n", r2))
