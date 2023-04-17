library(ggplot2)
library(gridExtra)

setwd("/Users/mfarman")

df <- read.table("TMRREs.txt", header = FALSE)

colnames(df) <- c('lineage', 'age')

ggplot(df, aes(x=age)) + geom_bar(color = "blue", fill = "blue") + 
  facet_grid(vars(lineage), scale = "free") + xlim(-5,10000) +
  ggtitle("Distribution of Time to Most Recent Recombination Event")+
  xlab("years before present") +
  ylab("# of interations") +
  ylim(-1, 201)
