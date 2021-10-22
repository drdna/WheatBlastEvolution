#Author: Mostafa Rahnama 
# Plotting  1) numbers of called sites on each chromosome; 2) # sites called for each population on each chromosome

library(readxl)
library(ggplot2)
library(tidyverse)
library(tidyverse)
library(rio)

sites <- read.csv("~/Google Drive/1.NATURE_GENETICS/FIGS5B_TOTAL_GENOME_CONTRIB/CP_called_sites_by_chromosome.csv", header = FALSE)
colnames(sites) <- c("Chr", "Total")

Calls <- read.csv("/Users/mostafa/Google Drive/1.NATURE_GENETICS/FIGS5B_TOTAL_GENOME_CONTRIB/CP_called_donors_by_chromosomes.csv")

#chromosome equivalents = sites called for pop/total sites
dm <- merge(Calls, sites, by = "Chr")
#sum all values of each pop regardless of group
dm$value <- dm$X.sites/dm$Total
dm <- subset(dm , group == "PoL-PoT")

level_order <- factor(dm$Chr, levels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7"))
dm$variable <- factor(dm$pop, levels = c("E", "Lu", "St","U1","X"))

colors <- c(E = "#A3A500", Lu= "#CC3333", St= "#F564E3", U1 = "#F8766D", X= "#619CFF")
pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/FIGS5B_TOTAL_GENOME_CONTRIB/Fig5C.", 
               ".pdf", sep = ''), 6, 4)  # should be ~ 11 wide x 8.5 tall 
print(
ggplot(dm, aes(x = level_order, y = value, fill = variable, group = value)) + 
  geom_bar(position = position_stack(), stat = "identity") + 
  labs(title = " ", y = "Chromosome equivalents\n", x= " ",  fill = "Population")+
  theme(#plot.title = element_text(hjust = 1.5),
        legend.title=element_text(face="bold", size=10), 
        legend.text=element_text(face="bold", size=9),
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom", 
        legend.justification="bottom",
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,2, "cm"),
        axis.title.x = element_text(face="bold", size = 13),  # , margin = margin(1, 0, 40, 0), vjust = -3
        axis.title.y = element_text(face="bold", size = 11),  
        axis.text.x=element_text(size=8, face="bold", colour = "black", angle=0),
        axis.text.y = element_text(face="bold", size=7, colour = "black"),
        #axis.ticks.y = element_blank(),    #element_line(size=1),
        axis.ticks.length = unit(0.1,"cm"),
        axis.ticks = element_line(size = 0.3),

        strip.text.y = element_text(size=10, face = "bold", margin = margin(0,1,0,1, "cm"))) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
  scale_fill_manual(values=colors)

)
dev.off()

################ Panel D
dm2 <- aggregate(value ~ variable , dm, mean)
dm2$variable <- factor(dm2$variable, levels = c("E", "Lu", "St","U1","X"))

#
Calls4 <- merge(Calls, sites, by = "Chr")
#sum all values of each pop regardless of group
Calls4$value <- Calls4$X.sites/Calls4$Total
Calls5 <- subset(Calls4, select = -c(Chr, X.sites, Total))
Calls6 <- aggregate(value ~ group+pop , Calls5, mean)

Calls7 <- spread(Calls6, pop, value)
pdf(file=paste("/Users/mostafa/Google Drive/1.NATURE_GENETICS/FIGS5B_TOTAL_GENOME_CONTRIB/Fig5D", 
               ".pdf", sep = ''), 6, 4)  # should be ~ 11 wide x 8.5 tall 
print(
  ggplot(dm2, aes(x = variable, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    geom_bar(data = subset(Calls6, group == "PoL-PoT"), mapping =aes(x= "PoL-PoT",y = value, fill = pop, group = value), position = position_stack(), stat = "identity") +
    geom_bar(data = subset(Calls6, group == "PoT"), mapping =aes(x= "PoT",y = value, fill = pop, group = value), position = position_stack(), stat = "identity") +
    geom_bar(data = subset(Calls6, group == "PoL"), mapping =aes(x= "PoL",y = value, fill = pop, group = value), position = position_stack(), stat = "identity") +
    scale_x_discrete(limits = c("E", "Lu", "St","U1","X", "PoL-PoT", "PoT", "PoL")) +
    labs(title = " ", y = "Genome equivalents", x= "Donor populations\n",  fill = "Population")+
    theme(#plot.title = element_text(hjust = 1.5),
      plot.margin = unit(c(0,0,0,0.1), "cm"),
      legend.title=element_text(face="bold", size=10), 
      legend.text=element_text(face="bold", size=9),
      legend.direction = "horizontal",
      legend.key = element_rect(fill = "white"),
      legend.key.size = unit(0.2, "cm"),
      legend.position = "bottom", 
      legend.justification="bottom",
      axis.title.x = element_text(face="bold", size = 11, vjust = -1),  # , margin = margin(1, 0, 40, 0), vjust = -3
      axis.title.y = element_text(face="bold", size = 11, vjust = 2),  
      axis.text.x=element_text(size=7, face="bold", colour = "black", angle=0),
      axis.text.y = element_text(face="bold", size=7, colour = "black"),
      #axis.ticks.y = element_blank(),    #element_line(size=1),
      axis.ticks.length = unit(0.1,"cm"),
      axis.ticks = element_line(size = 0.2),
      
      strip.text.y = element_text(size=10, face = "bold", margin = margin(0,0.1,0,0.1, "cm"))) + 
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
   + scale_fill_manual(values=colors)
  
)
dev.off()




