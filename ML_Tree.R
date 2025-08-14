# Load required libraries
library(ape)
library(phangorn)

# Read aligned FASTA (should be aligned already)
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# Convert alignment to phyDat for ML analysis
phy_data <- phyDat(alignment, type = "DNA")

# Build initial NJ tree as starting topology
dm <- dist.ml(phy_data)
start_tree <- NJ(dm)

# Estimate ML tree using GTR + Γ model
fit <- pml(start_tree, phy_data)
fit_opt <- optim.pml(fit, model = "GTR", optGamma = TRUE, optInv = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

# Bootstrap ML tree (100 replicates)
set.seed(123)
bs_ml <- bootstrap.pml(fit_opt, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

# Compute bootstrap support values
ml_tree <- fit_opt$tree
ml_tree$node.label <-  prop.clades(ml_tree, bs_ml)

# Rename tip labels
label_map <- c(
  "699985" = "A4ac", "760355" = "A4eo", "33881" = "A3bt", "22415" = "A3bb",
  "22385" = "A3", "34175" = "A3fa", "34422" = "A3u", "36733" = "A3ah",
  "35231" = "A3ad", "502153" = "A3ez", "363640" = "A3ey", "363655" = "A3fb",
  "22469" = "A3g"
)
ml_tree$tip.label <- sapply(ml_tree$tip.label, function(x) {
  if (!is.null(label_map[[x]])) label_map[[x]] else x
})

ntips <- Ntip(ml_tree)

# Plot the ML tree with bootstrap values ≥ 50
par(mar = c(6, 4, 4, 2))
plot(unroot(ml_tree), show.tip.label = TRUE, cex = 0.8,
     main = "A) Maximum Likelihood Tree of ITS2 Types with Bootstrap Support ≥50%", y.lim = c(-2, ntips))
nodelabels(text = ifelse(as.numeric(ml_tree$node.label) >= 50, ml_tree$node.label, ""),
           frame = "rect", cex = 1, font = 2)

plot(unroot(ml_tree), 
     show.tip.label = TRUE, 
     cex = 1.2, 
     font = 1,  # plain style
     family = "Times New Roman",  # or "Times New Roman"
     main = "A) Maximum Likelihood Tree of ITS2 Types with Bootstrap Support ≥50%",
     y.lim = c(-2, ntips))
# Extract bootstrap values and node positions
node_labels <- as.numeric(ml_tree$node.label)
label_positions <- which(node_labels >= 50)

# Add labels only for nodes with support ≥ 50, with background
nodelabels(text = node_labels[label_positions],
           node = label_positions + Ntip(ml_tree),  # shift to internal node numbers
           frame = "rect", cex = 1, font = 5)


# Custom scale bar
scale_length <- 0.07  # in tree branch length units
x_start <- usr[1] + 0.003
y_start <- usr[3] + 2.7

# Draw scale bar manually
segments(x_start, y_start, x_start + scale_length, y_start, lwd = 2)
text(x = x_start + scale_length / 2, y = y_start - 0.3,
     labels = "0.07 substitutions/site (model-based ML)", cex = 1, family = 'Times')
