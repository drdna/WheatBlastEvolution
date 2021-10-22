
##Generating plots for portion of chr length that has SNPs
## Authour: Mostafa Rahnama 

library(readxl)
library(ggplot2)
library(tidyverse)
library(rio)

# preparing datasets from the /CPApr042921/Out_3/B71.chr*.V2.complete_out.chunklengths.out in an excel file 
data_list <- import_list("~/Google Drive/1.NATURE_GENETICS/FigS7/B71.ALL_CHRs.V2.complete_out.chunklengths.out.xlsx", setclass = "tbl", rbind = TRUE)
colnames(data_list)[colnames(data_list) == '_file'] <- 'Chrs'
data_list$Chrs <- paste("Chr", data_list$Chrs, sep = "")
data_list$Recipient <- ifelse(data_list$Recipient ==  "pg1213-22", "Pg1213-22", data_list$Recipient)
data_list$Recipient <- ifelse(data_list$Recipient ==  "Br80", "BR80", data_list$Recipient)
data_list$Recipient <- ifelse(data_list$Recipient ==  "Br130", "BR130", data_list$Recipient)

strains <- read_excel("~/Google Drive/1.NATURE_GENETICS/Fig4/Final Figure/New_Genotypes4.xlsx", sheet = "Sheet1")
strains <- select(strains, c(2,3,15))
#strains$geneo <- paste(strains$Host, strains$genotype, sep = " ")

my_data1 <- subset(data_list, Recipient %in% strains$Recipient)
my_data2 <- merge(my_data1, strains, by.y = "Recipient")

# adding E1, E2 and E3 to E
my_data2 <- my_data2 %>% mutate("E" = E1+ E2+E3) 

#Lee+P+O+S (rename this last one X)
my_data3 <- my_data2 %>% mutate(X = Lee+P+O+S+P2)
my_data3 <- subset(my_data3, select= c("E","Lu", "St", "U1","X", "Chrs", "Geno"))   #"geneo", 

#write.csv(my_data3, "~/Google Drive/1.NATURE_GENETICS/FigS7/FigS7A_Chrs_Proportions.csv", row.names = FALSE) 

dm = reshape2::melt(my_data3, id.vars = c('Geno', 'Chrs'))

#level_order <- factor(dm$Geno, levels = strains$Geno)
level_order <- factor(dm$Geno, levels = c( "PoL1","PoT1","PoL2","PoL3","PoL4","PoL5","PoL6","PoL7","PoL8","PoL9","PoL10","PoL11",          "PoL13",  
                                          "PoT2","PoT3","PoT4","PoT5","PoT6","PoT7","PoT8","PoT9","PoT10","PoT11",
                                          "PoT12","PoT13","PoT14","PoT15","PoT16","PoT17","PoT18","PoT19","PoT20","PoT21","PoT22",
                                          "PoT24","PoT25","PoT26","PoT27","PoT28",
                                          "PoT30","PoT31","PoT32",          "PoT34") )


dm$variable <- factor(dm$variable, levels = c("E",  "Lu", "St","U1","X"))

colors <- c(E = "#A3A500", Lu= "#CC3333", St= "#F564E3", U1 = "#F8766D", X= "#619CFF")  #E3="#717200", 


          
pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/FigS7/Final_Figures/Final_Final_FigS7_A", 
               ".pdf", sep = ''), 6.5, 5)  # should be ~ 11 wide x 8.5 tall 
print(
ggplot(dm, aes(x = level_order, y = value*100, fill = variable, group = value)) + 
  geom_bar(position = "fill", stat = "identity") + 
  #scale_fill_discrete(breaks = c("C1","C2","'E1+E2'","E3","Ec", "Er","L2","Lee","Lu","O","P","S","St","U1","U2") ) +
  labs(title = " ", y = "% Core chromosomes from donor", x= "Haplotype\n",  fill = "Population")+
  scale_y_continuous(breaks = c(0, .50, 1.00),
  labels = function(x) paste0(round(as.numeric(x*100)), ""))+
  facet_grid(Chrs ~ ., scales = "free", space = "free")+
  theme(#plot.title = element_text(hjust = 1.5),
        legend.title=element_text(face="bold", size=5), 
        legend.text=element_text(face="bold", size=4.5),
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(0.15, "cm"),
        legend.position = "bottom", 
        legend.justification="bottom",
        #legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(0,0,0,2, "cm"),
        axis.title.x = element_text(face="bold", size = 6.5),  # , margin = margin(1, 0, 40, 0), vjust = -3
        axis.title.y = element_text(face="bold", size = 6.5),  
        axis.text.x=element_text(size=4, face="bold", colour = "black", angle=45, vjust=1, hjust = 1),
        axis.text.y = element_text(face="bold", size=3.5, colour = "black"),
        axis.ticks.y = element_blank(),    #element_line(size=1),
       # axis.ticks.length.x = unit(1,"cm"),
        axis.ticks.x = element_line(size = 0.15),

        strip.text.y = element_text(size=5, face = "bold", margin = margin(0,0.1,0,0.1, "cm"))) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
  scale_fill_manual(values=colors)

)
dev.off()

############################# Add up all chromosomes together 
dm2 <- aggregate(value ~ variable+Geno , dm, sum)
TABLES3 <- spread(dm2, key = variable, value = value)

write.csv(TABLES3, "~/Google Drive/1.NATURE_GENETICS/FigS7/FigS7B_Whole_Proportions_2.csv", row.names = FALSE) 
#level_order <- factor(dm2$Geno, levels = strains$Geno)
level_order <- factor(dm2$Geno, levels = c( "PoL1","PoT1","PoL2","PoL3","PoL4","PoL5","PoL6","PoL7","PoL8","PoL9","PoL10","PoL11",          "PoL13",  
                                           "PoT2","PoT3","PoT4","PoT5","PoT6","PoT7","PoT8","PoT9","PoT10","PoT11",
                                           "PoT12","PoT13","PoT14","PoT15","PoT16","PoT17","PoT18","PoT19","PoT20","PoT21","PoT22",
                                           "PoT24","PoT25","PoT26","PoT27","PoT28",
                                           "PoT30","PoT31","PoT32",          "PoT34") )

dm2$variable <- factor(dm2$variable, levels = c("E", "Lu", "St","U1","X"))

pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/FigS7/Final_Figures/Final_Final_FigS7_B", 
               ".pdf", sep = ''), 6.5, 2)  # should be ~ 11 wide x 8.5 tall 
print(
  ggplot(dm2, aes(x = level_order, y = value, fill = variable, group = value)) + 
    geom_bar(position = "fill",stat = "identity") + 
    #scale_fill_discrete(breaks = c("C1","C2","'E1+E2'","E3","Ec", "Er","L2","Lee","Lu","O","P","S","St","U1","U2") ) +
    labs(title = " ", y = "% Core genome from donor", x= "Haplotype\n",  fill = "Population")+
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                      labels = function(x) paste0(round(as.numeric(x*100)), ""))+
    #scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 10), labels = c("0", "",  "", "", "","1", "", "", "", "", "2", "", "")) +
    #facet_grid(Chrs ~ ., scales = "free", space = "free")+
    theme(#plot.title = element_text(hjust = 1.5),
      legend.title=element_text(face="bold", size=5), 
      legend.text=element_text(face="bold", size=4.5),
      legend.direction = "horizontal",
      legend.key = element_rect(fill = "white"),
      legend.key.size = unit(0.15, "cm"),
      legend.position = "bottom", 
      legend.justification="bottom",
      #legend.margin=margin(0,0,0,0),
      #legend.box.margin=margin(0,0,0,2, "cm"),
      axis.title.x = element_text(face="bold", size = 6.5),  # , margin = margin(1, 0, 40, 0), vjust = -3
      axis.title.y = element_text(face="bold", size = 6.5),  
      axis.text.x=element_text(size=4, face="bold", colour = "black", angle=45, vjust=1, hjust = 1),
      axis.text.y = element_text(face="bold", size=3.5, colour = "black"),
      axis.ticks.y = element_blank(),    #element_line(size=1),
      # axis.ticks.length.x = unit(1,"cm"),
      axis.ticks.x = element_line(size = 0.15),
      
      strip.text.y = element_text(size=5, face = "bold", margin = margin(0,0.1,0,0.1, "cm"))) + 
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
    scale_fill_manual(values=colors)
  
)


dev.off()






