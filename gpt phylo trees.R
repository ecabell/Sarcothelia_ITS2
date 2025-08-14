###############################
## 1.  LOAD PACKAGES
###############################
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(readr)

###############################
## 2.  READ THE ALIGNMENT
###############################
dna_aln <- read.FASTA("~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta")
dna_phydat <- phyDat(dna_aln, type = "DNA")

###############################
## 3.  STARTING TREE (NEIGHBOUR-JOINING)
###############################
d <- dist.hamming(dna_aln)
nj_start <- NJ(d)

###############################
## 4.  MAXIMUM PARSIMONY SEARCH
###############################
mp_tree <- pratchet(dna_phydat, start = nj_start, trace = 0)
mp_tree <- acctran(mp_tree, dna_phydat)

###############################
## 5.  BOOTSTRAP
###############################
set.seed(2025)
bs_reps <- 1000

bs_trees <- bootstrap.phyDat(
  dna_phydat,
  FUN   = pratchet,
  bs    = bs_reps,
  trace = 0
)

# Add bootstrap support
bs_vals <- round(100 * prop.clades(mp_tree, bs_trees) / bs_reps)
mp_tree$node.label <- bs_vals

# Optional: hide weakly supported nodes
mp_tree$node.label[bs_vals < 50] <- ""

###############################
## 6.  Export Tree
###############################
write.tree(mp_tree, file = "MP_bootstrap_final.nwk")
write.nexus(mp_tree, file = "MP_bootstrap_final.nex")

###############################
## 7.  Load metadata
###############################
meta <- read_csv("~/Documents/Zoox PCR/meta.csv", col_types = "cccc")
mp_tree$tip.label <- as.character(mp_tree$tip.label)

###############################
## 8.  Plot MP tree with ggtree (various layouts)
###############################

layouts <- c("fan", "roundrect", "circular", "unrooted", "daylight")
for (layout in layouts) {
  p <- ggtree(mp_tree, layout = layout) %<+% meta +
    geom_tippoint(aes(shape = Location, colour = Morph), size = 4) +
    geom_tiplab(aes(label = ITS2_Type), size = 5, hjust = -0.2) +
    geom_label2(aes(
      subset = !isTip & !is.na(label),
      label = label
    ),
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    label.r = unit(0.2, "lines"),
    fill = "grey95",
    color = 'grey40',
    size = 3.5
    ) +
    scale_colour_manual(values = c(
      "Blue"        = "blue",
      "Brown"       = "#D55E00",
      "Blue, Brown" = "#CC79A7"
    )) +
    scale_shape_manual(values = c(
      "Kauai"                   = 15,
      "Oahu"                    = 16,
      "Big Island"              = 18,
      "Big Island, Oahu"        = 17,
      "Big Island, Kauai"       = 25,
      "Big Island, Oahu, Kauai" = 8
    )) +
    ggtitle(paste("Maximum-Parsimony Tree with Bootstrap –", layout)) +
    theme_tree2() + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.key.size = unit(1, "cm"))
  
  print(p)
  ggsave(paste0("MP_tree_", layout, ".pdf"), p, width = 10, height = 8)
}



########NJ

###############################
## NJ TREE FROM GENETIC DISTANCES
###############################

# Compute model-based distances (substitutions/site)
genetic_dist <- dist.ml(dna_phydat)

# Build NJ tree
nj_tree <- NJ(genetic_dist)

# Bootstrap NJ tree
set.seed(2025)
bs_nj_trees <- bootstrap.phyDat(
  dna_phydat,
  FUN = function(x) NJ(dist.ml(x)),
  bs = 1000
)

# Bootstrap values
nj_bs_vals <- round(100 * prop.clades(nj_tree, bs_nj_trees) / 1000)
nj_tree$node.label <- nj_bs_vals
nj_tree$node.label[nj_bs_vals < 50] <- ""

# Export tree
write.tree(nj_tree, file = "NJ_bootstrap_final.nwk")
write.nexus(nj_tree, file = "NJ_bootstrap_final.nex")

for (layout in layouts) {
  p <- ggtree(nj_tree, layout = layout) %<+% meta +  # changed from mp_tree
    geom_tippoint(aes(shape = Location, colour = Morph), size = 4) +
    geom_tiplab(aes(label = ITS2_Type), size = 5, hjust = -0.2) +
    geom_label2(aes(
      subset = !isTip & !is.na(label),
      label = label
    ),
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    label.r = unit(0.2, "lines"),
    fill = "grey95",
    color = 'grey40',
    size = 3.5
    ) +
    scale_colour_manual(values = c(
      "Blue"        = "blue",
      "Brown"       = "#D55E00",
      "Blue, Brown" = "#CC79A7"
    )) +
    scale_shape_manual(values = c(
      "Kauai"                   = 15,
      "Oahu"                    = 16,
      "Big Island"              = 18,
      "Big Island, Oahu"        = 17,
      "Big Island, Kauai"       = 25,
      "Big Island, Oahu, Kauai" = 8
    )) +
    scale_x_continuous(name = "Genetic Distance (substitutions/site)") +  # adds correct scale
    ggtitle(paste("NJ Tree with Bootstrap –", layout)) +
    theme_tree2() +
    theme(
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 15),
      legend.key.size  = unit(1, "cm")
    )
  
  print(p)
  ggsave(paste0("NJ_tree_", layout, ".pdf"), p, width = 10, height = 8)
}


#######ML tree 

###############################
## ML TREE USING GTR MODEL
###############################

# Initial NJ tree as starting topology
d_ml <- dist.ml(dna_phydat)
nj_start <- NJ(d_ml)

# Initial ML fit
fit <- pml(nj_start, data = dna_phydat)

# Optimize tree and model parameters
fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

ml_tree <- fit_opt$tree  # Extract final ML tree

###############################
## ML BOOTSTRAP
###############################
set.seed(2025)
bs_ml <- bootstrap.pml(fit_opt, bs = 1000, optNni = TRUE, multicore = FALSE)

ml_tree$node.label <- round(100 * prop.clades(ml_tree, bs_ml) / 1000)
ml_tree$node.label[ml_tree$node.label < 50] <- ""



for (layout in layouts) {
  p <- ggtree(ml_tree, layout = layout) %<+% meta +
    geom_tippoint(aes(shape = Location, colour = Morph), size = 4) +
    geom_tiplab(aes(label = ITS2_Type), size = 5, hjust = -0.2) +
    geom_label2(aes(
      subset = !isTip & !is.na(label),
      label = label
    ),
    label.padding = unit(0.2, "lines"),
    label.size = 0.4,
    label.r = unit(0.2, "lines"),
    fill = "grey95",
    color = 'grey40',
    size = 3.5
    ) +
    scale_colour_manual(values = c(
      "Blue"        = "blue",
      "Brown"       = "#D55E00",
      "Blue, Brown" = "#CC79A7"
    )) +
    scale_shape_manual(values = c(
      "Kauai"                   = 15,
      "Oahu"                    = 16,
      "Big Island"              = 18,
      "Big Island, Oahu"        = 17,
      "Big Island, Kauai"       = 25,
      "Big Island, Oahu, Kauai" = 8
    )) +
    scale_x_continuous(name = "Genetic Distance (substitutions/site)") +
    ggtitle(paste("ML Tree with Bootstrap –", layout)) +
    theme_tree2() +
    theme(
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 15),
      legend.key.size  = unit(1, "cm")
    )
  
  print(p)
  ggsave(paste0("ML_tree_", layout, ".pdf"), p, width = 10, height = 8)
}

fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))
logLik(fit_opt)
AIC(fit_opt)
BIC(fit_opt)



- **Compare trees** (e.g., NJ vs ML with same model and data):
  ```r
logLik(fit_opt)  # ML tree
logLik(pml(nj_tree, data = dna_phydat))  # NJ tree
AIC(pml(nj_tree, data = dna_phydat))
BIC(pml(nj_tree, data = dna_phydat))
