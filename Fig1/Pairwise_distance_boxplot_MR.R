library(ggplot2)
library(reshape2)
library("PopGenome")

#******************** Can you modify this R-CODE to include horizontally jittered dots for each value superimposed on the box plots.


my_df <- read.table("/Users/mostafa/Google Drive/1.NATURE_GENETICS/Fig1D/boxplot-distances_CC2.txt", header=FALSE)

colnames(my_df) <- c("Q", "S", "distance", "population")
#my_df <- melt(my_df)
#my_df$distance <- as.numeric(my_df$distance)
my_df$population <- as.factor(my_df$population)
#means2 <- aggregate( ~ population, my_df, function(x) c(mean = mean(x), sd = sd(x)))
means2 <- my_df %>% group_by(population) %>% 
                              summarise(Numbers = n(), Mean = mean(distance/1000), SD = sd(distance/1000))
write.csv(means2,"~/Google Drive/1.NATURE_GENETICS/1.SUBMISSION/SupplTables/Distance_table_Fig1D", row.names = FALSE)

#my_df <- subset(my_df, population != "Er")
#add a dummy number for Er in the range of 
#my_df[nrow(my_df)+1, ] = list(Q="EtKY19-1",S="G17",distance=as.numeric(386),population="Er")


#my_df$colr <- factor(my_df$population, levels = unique(my_df$population))

#colors <- c(Er="#CC6666", L2="#666600", Lu= "#CC3333", 
#           Ec="#669999",  C2= "#FFA500", E3="#717200",
#          U2 = "#D2B48C", C1="#9370DB", Lee="#C0C0C0", O="#00CCCC",  
#         S= "#619CFF", 'E1+E2' = "#A3A500", U1 = "#F8766D", P = "#00CC66", St = "#F564E3")

colors <- c(C1="#FFFF00", E1= "#A3A500", E2="#717200",Ec="#669999", Er="#CC6666", L="#9370db", Lee="#C0C0C0", Lu= "#CC3333",O="#00CCCC",
            P="#00CC66", S= "#619CFF",St = "#F564E3", T= "#0000FF" , U1 = "#F8766D", U3 = "#D2B48C", U4= "#FF8C00")
#new pops: E1 E2,L, T,U3, U4


#pdf("Desktop/2.SCIENCE2020/FigS1.pdf", 11.5, 6.5)
pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/Fig1D/Fig1D_5", 
               ".pdf", sep = ''),  6.5, 3.5, useDingbats=FALSE)
print(
ggplot(my_df, aes(x = population, y = distance/1000000)) +
  geom_jitter(aes(color= factor(population)), alpha = 0.5) + #
  
  geom_boxplot(lwd=0.5, outlier.shape = NA) + 
 # geom_point(position=position_jitterdodge(jitter.width=2, dodge.width = 0), 
  #           pch=21, aes(fill=factor(population), ), show.legend = F)+
 # geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
  #           aes(color=factor(population)), show.legend = F) +
  #
  #geom_text(data = means, aes(label = distance/1000000)) +
  scale_fill_manual(name="", values=colors)+
  scale_color_manual(name="", values=colors)+
  
  theme_bw() + 
  xlab("Population") +
  ylab("Nucleotide Diversity") +
  theme(#plot.title = element_text(hjust = 1.5),
    legend.title=element_text(face="bold", size=5), 
    legend.text=element_text(face="bold", size=8),
    legend.direction = "vertical",
    legend.key = element_rect(fill = "white"),
   # legend.key.size = unit(1, "cm"),
    legend.position = "right", 
    #legend.justification="right",
    #legend.margin=margin(0,0,0,0),
    #legend.box.margin=margin(0,0,0,2, "cm"),
    axis.title = element_text(face="bold", size = 13),  # , margin = margin(1, 0, 40, 0), vjust = -3
    #axis.title.y = element_text(face="bold", size = 10),  
    axis.text=element_text(size=8, face="bold", colour = "black", angle=0),
   panel.grid.major.y = element_blank(),
   panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
   panel.grid.minor.x = element_blank()
    #axis.text.y = element_text(face="bold", size=3.5, colour = "black"),
   # axis.ticks.y = element_blank(),    #element_line(size=1),
    # axis.ticks.length.x = unit(1,"cm"),
    #axis.ticks.x = element_line(size = 0.15)
  )
#+  guides(colour=guide_legend(nrow=2)) 
+ guides(color = guide_legend(override.aes = list(size=5)))

)
dev.off()



######## Test if there is any sig differences between species' distances, between the different Pops
# using Kruskal-Wallis test because the requirement of one-way ANOVA don't meet here (because data are not normal)

#significance test
kr <- kruskal.test(distance ~ population, my_df) #Significant differences as determined by Kruskal-Wallis rank sum test
mw <- pairwise.wilcox.test(my_df$distance, my_df$population)  ## Pairwise Wilcoxon (or "Mann-Whitney") rank sum tests between all factor level combinations

#using "pairwise significance grouping labels" for finding similar letters
#matrix showing connections between levels  https://www.r-bloggers.com/2014/05/automated-determination-of-distribution-groupings-a-stackoverflow-collaboration/
n = 16  # number of group levels to tes

# Create matrix showing factor levels that should be grouped
g <- as.matrix(mw$p.value > 0.05) # TRUE signifies that pairs of not significantly different at the p < 0.05 level
g <- cbind(rbind(NA, g), NA) # make square
g <- replace(g, is.na(g), FALSE) # replace NAs with FALSE
g <- g + t(g) # not necessary, but make matrix symmetric
diag(g) <- 1 # diagonal equals 1
rownames(g) <- 1:n # change row names
colnames(g) <- 1:n # change column names
g # resulting matrix

library(igraph)
# Re-arrange data into an "edge list" for use in igraph (i.e. which groups are "connected") - Solution from "David Eisenstat" ()
same <- which(g==1)
g2 <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
g2 <- g2[order(g2[[1]]),] # Get rid of loops and ensure right naming of vertices
g3 <- simplify(graph.data.frame(g2,directed = FALSE))
get.data.frame(g3) # view connections

# Plot igraph
png("/Users/mostafa/Google Drive/1.NATURE_GENETICS/Fig1D/igraph_level_groupings_0.5.png", width=5, height=5, units="in", res=400, type="cairo")
par(mar=c(3,1,1,1))
V(g3)$color <- 8
V(g3)$label.color <- 1
V(g3)$size <- 20
plot(g3) # plot all nodes are connections
box()
mtext("Linked levels are not significantly different \n(Mann-Whitney)", side=1, line=1)
dev.off()

# Calcuate the maximal cliques - these are groupings where every node is connected to all others
cliq <- maximal.cliques(g3) # Solution from "majom" ()

# Reorder by level order - Solution from "MrFlick" ()
ml<-max(sapply(cliq, length))
#reord <- do.call(order, data.frame(  do.call(rbind, 
  #        lapply(cliq, function(x) c(sort(x), rep.int(0, ml-length(x))))))

sf <- function(x) sort(x)[seq_len(ml)]
cliq2 <- lapply(cliq, as.numeric)
reord <- do.call(order, as.data.frame(do.call(rbind, lapply(cliq2, sf))))

cliq <- cliq[reord]
cliq
