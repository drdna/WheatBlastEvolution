## Plots trees and colors taxa according to lineage affinities


library(ape)
library(phangorn)
library(phytools)
library(RColorBrewer)
library(tidyverse)
library(ggtree)

#read in data
Tree1 <- read.tree("PoTPoL_tree2.nwk")

# assign the groups
groupInfo <- split(Tree1$tip.label, gsub(".*_", "", Tree1$tip.label))

Tree1 <- groupOTU(Tree1, groupInfo)

Tree1$tip.label <- gsub("_.*", "", Tree1$tip.label)

# create the palette
n <- 13
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#colorVector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#colorVector = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#myColors <- sample(colorVector, n)

myColors <- c(Poa = "black", Triticum= "#0000FF",  Cenchrus="#00FF00", Elionurus = "red",
            Eleusine= "#A3A500", 
            Lolium = "#a64dff", Digitaria="#666600", Hordeum= "darkgray", Bromus = "#008000",
            Melinis= "#00BFC4", Avena= "#FF00FF", Urochloa= "#CC0000", Echinochloa = "orange")

# plot and save the tree
options(ignore.negative.edge=TRUE)

Tree2 <- drop.tip(Tree1, "BdBar", subtree=FALSE)

#pdf("PoT_PoL_tree.pdf", 8.5, 11)
ggtree(Tree1, 
       layout = "circular", 
       branch.length='none', 
       aes(color = group, label = gsub("_.*", "", label)))+
 geom_tiplab(size=4, offset = 1) +
  scale_color_manual(values = myColors) +
  guides(col = guide_legend(nrow = 13, width = 5, override.aes = aes(size = 2, label = ""))) +
  labs(color = "Host Genera")
#dev.off()
#ggsave("tree.png", dpi = 600, width =10, height =10)

geom_tippoint(aes(color=location), size=3, alpha=.75) +
  scale_color_brewer("location", palette="Spectral")




