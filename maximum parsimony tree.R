library(ape)
library(phangorn)


###Ultimately opted out of using MP tree in study because can't show genetic distance
###my goal is to show genetic distance


# Read aligned FASTA (should be aligned already)
alignment <- read.dna("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta", format = "fasta")

# Convert to phyDat object for MP analysis
phy_data <- phyDat(alignment, type = "DNA")

# Compute distance matrix
dm <- dist.ml(phy_data)

# Generate Neighbor-Joining tree as a starting point
nj_tree <- NJ(dm)

# Optimize the NJ tree using parsimony
mp_tree <- optim.parsimony(nj_tree, phy_data)

# --- Bootstrap MP tree (correct method) ---
set.seed(123)  # For reproducibility
bs_mp <- bootstrap.phyDat(phy_data, FUN = function(x) {
  dm_bs <- dist.ml(x)
  tree_start <- NJ(dm_bs)
  optim.parsimony(tree_start, x)
}, bs = 100)

bs_vals <- prop.clades(mp_tree, bs_mp)

# Assign these to the internal node labels of your tree
mp_tree$node.label <- bs_vals

label_map <- c(
  "699985" = "A4ac",
  "760355" = "A4eo",
  "33881" = "A3bt",
  "22415" = "A3bb",
  "22385" = "A3",
  "34175" = "A3fa",
  "34422" = "A3u",
  "36733" = "A3ah",
  "35231" = "A3ad",
  "502153" = "A3ez",
  "363640" = "A3ey",
  "363655" = "A3fb",
  "22469" = "A3g"
  # Add all mappings here
)


# Rename tip labels at the end
mp_tree$tip.label <- sapply(mp_tree$tip.label, function(x) {
  if (!is.null(label_map[[x]])) {
    label_map[[x]]
  } else {
    x
  }
})
# Plot and add node labels manually
plot(mp_tree, main = "Maximum Parsimony Tree with Bootstrap (%)")
nodelabels(bs_vals, frame = "none", cex = 0.7)
#edgelabels()


# Add bootstrap values to the tree
plotBS(mp_tree, bs_mp, p = 0, main = "Maximum Parsimony Tree with Bootstrap ≥ 50")

# Extract bootstrap values (majority rule consensus)
cons_tree <- consensus(bs_mp, p = 0.5)
bs_vals <- prop.clades(mp_tree, bs_mp)

# Plot and add node labels manually
plot(mp_tree, main = "Maximum Parsimony Tree with Bootstrap (%)")
nodelabels(bs_vals, frame = "none", cex = 0.7)
#edgelabels()

plotBS(mp_tree, bs_mp, p = 0, main = "All Bootstrap Values")
plotBS(mp_tree, bs_mp, p = 50, main = "Only Bootstrap ≥ 50")
plotBS(mp_tree, bs_mp, p = 90, main = "Only Bootstrap ≥ 90")


# Plot with node labels as bootstrap proportions
plot(cons_tree, main = "Majority-Rule Bootstrap Consensus Tree")
#nodelabels()


write.tree(mp_tree, file = "mp_tree_final.nwk")
pdf("mp_tree_final_plot.pdf", width = 7, height = 10)
dev.off()


label_map <- c(
  "699985" = "A4ac",
  "760355" = "A4eo",
  "33881" = "A3bt",
  "22415" = "A3bb",
  "22385" = "A3",
  "34175" = "A3fa",
  "34422" = "A3u",
  "36733" = "A3ah",
  "35231" = "A3ad",
  "502153" = "A3ez",
  "363640" = "A3ey",
  "363655" = "A3fb",
  "22469" = "A3g"
  # Add all mappings here
)


# Rename tip labels at the end
mp_tree$tip.label <- sapply(mp_tree$tip.label, function(x) {
  if (!is.null(label_map[[x]])) {
    label_map[[x]]
  } else {
    x
  }
})

