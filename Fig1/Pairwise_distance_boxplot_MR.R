library(ggplot2)
library(reshape2)
library(dplyr)

#******************** Can you modify this R-CODE to include horizontally jittered dots for each value superimposed on the box plots.


my_df <- read.table("~/Newest_boxplot.txt", header=FALSE)

colnames(my_df) <- c("Q", "S", "distance", "population")
my_df$population <- as.factor(my_df$population)

my_df <- my_df[my_df$Q != "87-120", ]
my_df <- my_df[my_df$S != "87-120", ]

#means2 <- aggregate( ~ population, my_df, function(x) c(mean = mean(x), sd = sd(x)))
means2 <- my_df %>% group_by(population) %>% summarize(Numbers = n(), Mean = mean(distance/1000), SD = sd(distance/1000))

colors <- c(PoC1="#FFFF00", PoE1= "#A3A500", PoE2="#717200",PoEc="#669999", PoEr1="#CC6666",
            PoL1="#9370db", PoLe="#C0C0C0", PoLu= "#CC3333",PoO="#00CCCC",
            PoP="#00CC66", PoS= "#619CFF",PoSt = "#F564E3", PoT= "#0000FF" , PoU1 = "#F8766D", PoU3 = "#D2B48C", PoU4= "#FF8C00")


pdf(file="~/Fig1D_5.pdf", 8, 5, useDingbats=FALSE)

print(ggplot(my_df, aes(x = population, y = distance/1000000)) +
  geom_boxplot(lwd=0.5, outlier.shape = NA) +  scale_fill_manual(name="", values=colors)+
  geom_jitter(aes(color= factor(population)), size = 4, alpha = 0.1, width = 0.2) +
  scale_color_manual(name="", values=colors)+
  theme_bw() + 
  ggtitle("Within-Lineage Nucleotide Diversity") +
  xlab("\nLineage") +
  ylab("Nucleotide Diversity (\u03c0)\n") +
  theme(#plot.title = element_text(hjust = 1.5),
  legend.title=element_text(face="bold", size=5), 
  legend.text=element_text(face="bold", size=8),
  legend.direction = "vertical",
  legend.key = element_rect(fill = "white"),
  legend.position = "none", 
  axis.title = element_text(face="bold", size = 14),
  axis.text=element_text(size=10, face="bold", colour = "black", angle=-45),
  axis.text.y = element_text(face="bold", size=10, colour = "black", angle = 0),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
  )
)

dev.off()
