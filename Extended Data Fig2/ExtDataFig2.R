# Description -------------------------------------------------------------

## Takes two trees and rotates branches to generate best alignment

## connects strain IDs with lines

library(ape)
library(phangorn)
library(phytools)

# Provide the objects

Tree1 <- read.tree("ch7bac7noGaps.tre")

Tree1

Tree2 <- read.tree("MPG1noGaps.tre")

cophylo <- cophylo(Tree1, Tree2, assoc=NULL, rotate=TRUE)

pdf("cophyloPlot.pdf", 24, 44)

plot(cophylo)

dev.off()

