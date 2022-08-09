# Root-to-tip regression

library(ggplot2)
library(ggpubr)

df <- as.data.frame(read.table("Root-to-tip2_clades.txt", header = T))

rownames(df) <- sub('_.+$', '', df$tip)

myColors <-(c("purple", "blue"))


p <- ggplot(df) + geom_point(aes(x=date, y=distance, col=clade, label = rownames(df)), size = 4) + geom_smooth(aes(x=date, y=distance, col=clade), method="lm", se=FALSE, fullrange=TRUE, show.legend = FALSE) +
  stat_regline_equation(label.x=1932, aes(x=date, y=distance, col=clade, label = ..rr.label..), size = 6, show.legend = FALSE) +
  ggtitle("Root to Tip Distances") +
  scale_color_manual(name = "Lineage", values = myColors, labels = c("Lolium", "Triticum")) +
  scale_x_continuous(expand = c(0, 0), limits=c(1930, 2020), breaks = c(1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 0.06)) + 
  geom_text() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x=element_text(size=12, vjust=-2),
        axis.text.y=element_text(size=12, hjust=2),
        axis.title.x=element_text(size = 14, face = "bold", vjust = -2),
        axis.title.y=element_text(size = 14, face = "bold", vjust = 2),
        plot.margin = margin(2,2,2,2, "cm")) +
  geom_segment(aes(x=1930,y=0,xend=2020,yend=0), color = "black") +
  geom_segment(aes(x=1930,y=0,xend=1930,yend=0.06), color = "black")

p + geom_smooth(aes(x=date, y=distance), color = "black", method="lm", se=FALSE, fullrange=TRUE, show.legend = FALSE) +
  stat_regline_equation(label.x=1932, label.y=0.052, aes(x=date, y=distance, label = ..rr.label..), size = 6, show.legend = FALSE)
  print (p)
  