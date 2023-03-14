# R
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ape)
library(scales)
library(patchwork)

# Load ML tree
t=read.tree("~/Chr1Chr2Chr5_undated.fasta.raxml.bestTree")

# Compute pairwise cophenetic / patristic distances and select distances to ATCC64557 - the oldest sample
distances <- cophenetic(t)
colnames(distances) = gsub("_.*", "", colnames(distances))

Dist2ATCC <- distances[colnames(distances) == 'ATCC64557', ]

# Load the collection dates and match them with the distances
dateLookup <- read.table('~/WB_dates.txt', header = FALSE)
dates <- c()
lineages <- c()

for(n in names(Dist2ATCC)){
  dates <- c(dates, dateLookup[dateLookup[,1] == n, 2])
  lineages <- c(lineages, dateLookup[dateLookup[,1] == n, 3])
}

df <- data.frame(Year = dates, Dist2ATCC = Dist2ATCC, lineages = lineages)
df <- df[rownames(df) != "ATCC64557", ]  # remove self comparison

# plot the tree

p <- ggscatter(df, x="Year", y="Dist2ATCC", size = 3, color = "lineages", add = "reg.line", 
               add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE, fullrange = TRUE) +
  scale_color_manual(name = "Lineage", values = c("purple", "blue"), labels = c(expression(italic("Lolium")), expression(italic("Triticum")))) +
  scale_x_continuous(name = "Sampling Date", expand = c(0, 0), limits=c(1980, 2020), breaks = c(1980, 1990, 2000, 2010, 2020)) +
  scale_y_continuous(name = "Patristic Distance to ATCC64557",expand = c(0, 0), limits=c(0, 0.08)) + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x=element_text(size=12, vjust=-2),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size = 12, face = "bold", vjust = -3),
        axis.title.y=element_text(size = 12, face = "bold", vjust = 2),
        plot.margin = margin(0,0,0,0, "cm"))

p1 <- p + stat_cor(method = "pearson", label.x = 1982, label.y = 0.075)


#plot(m, main = "Phylogenetic Signal for WB/GLS population)", xlim = c(1960, 2020),
#     xlab = "Collection Year", ylab = "Patristic Distance to ATCC64557")
#legend('topright', paste("Pearson\'s r =", round(cor(m$Coll_Year, m$Patr_dist_to_ATCC), 2)), bty = 'n')
#abline(lm(m$Patr_dist_to_ATCC ~ m$Coll_Year), lty = 2)
#cor.test(m$Coll_Year, m$Patr_dist_to_ATCC)

set.seed(123)
Resampled <- c()
Randomized <- c()

for(i in 1:1000) {
  
  resampledf <- df[sample(nrow(df), replace = TRUE), ]
  resampleCorr <- cor(resampledf$Year, resampledf$Dist2ATCC)
  Resampled <- c(Resampled, resampleCorr, method = "pearson")

  randomdf <- cbind(sample(df$Year, replace = FALSE), df$Dist2ATCC)
  randomCorr <- cor(randomdf[,1], randomdf[,2], method = "pearson")
  Randomized <- c(Randomized, randomCorr)
  
}

# melt data frame for plotting
Perms2 <- melt(cbind(Resampled, Randomized), id.vars = colnames) 

p2 <- ggplot(Perms2) + geom_boxplot(aes(x=Var2, y=value), color = "black") +
  xlab("Sampling Strategy") +
  ylab("R value") +
  theme(plot.title = element_text(size = 18, face = "bold", vjust = -2),
        axis.text.x=element_text(size=12, vjust = -2),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size = 12, face = "bold", vjust = -3),
        axis.title.y=element_text(size = 12, face = "bold", vjust = 2),
        plot.margin = margin(0,0,0,1, "cm"))

pdf("~/Fig5AB.pdf", 8.5, 4)
p1 + p2 + plot_layout(widths = c(5, 2.5)) + 
     plot_annotation(title = "Phylogenetic signal",
                     theme = theme(plot.title = element_text(hjust = 0.5, size = 18)))
dev.off()
#ggsave("Fig5AB.tiff", dpi = 600, bg = "white")
  