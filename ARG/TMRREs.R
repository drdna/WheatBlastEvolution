library(ggplot2)
library(gridExtra)

setwd("/Users/mfarman")

df <- read.table("TMRREsReal.txt", header = FALSE)

colnames(df) <- c('lineage', 'age')

ggplot(df, aes(x=age)) + geom_bar(color = "blue", fill = "blue") + 
  facet_grid(vars(lineage), scale = "free") + xlim(-1,1000)
