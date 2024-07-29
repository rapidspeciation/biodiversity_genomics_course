#exploring multiple trees
#
library(ggtree)
library(phytools)
library(dplyr)
library(ggdensitree)

TREE_FILE="concord.cf.tree"

#read in the concatenated tree
tree <- read.tree(TREE_FILE)

#create the plot
concat_tree <-
  ggtree(tree, ladderize = T) + 
  geom_tiplab() + 
  geom_treescale()

#take a look
concat_tree

#we can add some improvement
concat_tree <- 
  ggtree(tree, ladderize = T) + 
  geom_tiplab(align = T) + 
  geom_treescale() + 
  xlim(c(0,0.4)) +		#adjust the limits of the x axis so that we can see the tiplabels
  geom_nodelab(hjust = 1.1, vjust = -0.6) + #plots the node information (bootstrap and concordance factor)
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue") #add the branch length

concat_tree

#we can save the tree if we want
ggsave(plot=concat_tree, filename = "plots/concat_tree.png",
       height = 8,
       width = 6)


##Rerooting the tree
#check tiplabel or nodelabel index 
tree$tip.label

#reroot with phytools with outgroup
ggtree(phytools::reroot(tree,  node.number = 1), ladderize = T) + 
  geom_tiplab(align = T) + 
  geom_nodelab() +
  scale_x_continuous(limits = c(0,0.4))+
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

#specify the position of the root, 
#note that this position is arbitrary, we do not have the information unless we use one more outgroup
concat_tree_rerooted <- 
  ggtree(phytools::reroot(tree,  node.number = 1, position = 0.1), ladderize = T) + 
  geom_tiplab(align = T) + 
  geom_nodelab() +
  scale_x_continuous(limits = c(0,0.4))+
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

concat_tree_rerooted

ggsave(plot=concat_tree_rerooted, filename = "plots/concat_tree_rerooted.png",
       height = 8,
       width = 6)


###read in the file with multiple trees
trees <- read.tree("loci.treefile")

#create the plot
gene_trees <-
  ggdensitree(layout = "slanted", trees, color="lightblue") +
  geom_tiplab() +
  geom_treescale()

gene_trees

#set branch.length="none", this will render a cladogram
trees_slanted <- 
  ggdensitree(layout = "slanted", trees, branch.length="none", color="blue", alpha = 0.3) +
  geom_tiplab() +
  geom_treescale()

trees_slanted

#save the tree if you want
ggsave(plot = trees_slanted,
       device = "png",
       filename = "trees_slanted.png")

#plot the concatenated tree on top of the densitree
##change the concatenated tree to cladogram by setting branch.length "none" and layout "slanted"
concat_tree <- 
  ggtree(tree, ladderize = T, layout = "slanted", branch.length = "none") + 
  geom_tiplab(align = T) + 
  geom_treescale() + 
  #xlim(c(0,0.4)) +       
  geom_nodelab(hjust = 1.1, vjust = -0.6) +
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

#extract the data information to add to the densitree plot
d1 <- trees_slanted$data
d2 <- concat_tree$data
#move the start point of the root
d2$x <- d2$x - max(d2$x)
#check that they are the same
d1$x
d2$x

#Plot the trees
trees_slanted_conc <- 
trees_slanted +
  geom_tree(data = d2, layout = "slanted") +
  geom_nodelab(data = d2, hjust = 1.1, vjust = -0.6) +
  geom_text2(data = d2, aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "orange4")

ggsave(plot = trees_slanted_conc,
       device = "png",
       filename = "trees_slanted_conc.png")
