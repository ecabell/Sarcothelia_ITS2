#source maximum parsimony tree
library(ggtree)
library(tidyverse)

# Load metadata
meta <- read.csv("~/Desktop/meta.csv", stringsAsFactors = FALSE)

library(ggtree)
library(tidyverse)


ggplot(mp_tree) + geom_tree() + theme_tree() 
ggtree(mp_tree, color="firebrick", size=2, linetype="dotted")
ggtree(mp_tree, branch.length="none") + geom_tiplab()

ggtree(mp_tree)
ggtree(mp_tree, layout="roundrect") + geom_tiplab()
ggtree(mp_tree, layout="slanted") + geom_tiplab()
ggtree(mp_tree, layout="ellipse") + geom_tiplab()
ggtree(mp_tree, layout="circular") + geom_tiplab()
ggtree(mp_tree, layout="fan", open.angle=100) + geom_tiplab()
ggtree(mp_tree, layout="equal_angle") + geom_tiplab()
ggtree(mp_tree, layout="daylight") + geom_tiplab()
ggtree(mp_tree, branch.length='none') + geom_tiplab()
ggtree(mp_tree, layout="ellipse", branch.length="none") + geom_tiplab()
ggtree(mp_tree, branch.length='none', layout='circular') + geom_tiplab()
ggtree(mp_tree, layout="daylight", branch.length = 'none') + geom_tiplab()


ggtree(mp_tree) + layout_dendrogram()
ggplotify::as.ggplot(ggtree(mp_tree), angle=-30, scale=.9)
ggtree(mp_tree, layout='slanted') + coord_flip()
ggtree(mp_tree, layout='slanted', branch.length='none') + layout_dendrogram()
ggtree(mp_tree, layout='circular') + xlim(-10, NA) + geom_tiplab()
ggtree(mp_tree) + layout_inward_circular() + geom_tiplab()
ggtree(mp_tree) + layout_inward_circular(xlim=15)

ggtree(tree) + 
  geom_point(aes(shape=isTip, color=isTip), size=3)

p <- ggtree(mp_tree) + 
  geom_nodepoint(color="#b5e521", alpha=1/4, size=10) 
p + geom_tippoint(color="#FDAC4F", shape=8, size=5)
p + geom_tiplab(size=3, color="purple")

ggtree(mp_tree, layout="circular") + geom_tiplab(aes(angle=angle), color='black')

ggtree(mp_tree) + geom_tiplab(as_ylab=TRUE, color='firebrick')



