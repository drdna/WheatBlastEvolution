#Author: Mostafa Rahnama
# Plotting copyprobsperlocus.out data from CP based on chromosomes and strains of interest 

library(reshape2)
library(scales)
library(ggplot2)
library(grid)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(readxl)
library(readr)
library(tidyr)
'%notin%' <- Negate('%in%')

dt <- read.csv("~/Google Drive/1.NATURE_GENETICS/FigS9/2PoE1_example_on_Chr1.csv", header=TRUE)
dt <- subset(dt, select = -c(main))
dt1 <- subset(dt, pop %in% c("E1", "E2", "E3", "T") )

dt2 <- gather(dt1, sites, value, X4738637:X4874974)
dt2$sites <- gsub("X", "", as.character(dt2$sites))


### reading grouping data
#PoE2 <- read.table("~/Google Drive/1.NATURE_GENETICS/Fig XX/2PoE_groups")
#PoE2$groups <- 1:nrow(PoE2)
PoE2 <- read_excel("~/Google Drive/1.NATURE_GENETICS/FigS9/2PoE_groups.xlsx")
PoE2_1 <- gather(PoE2, coln, strain, V1:V56)
PoE2_2 <- PoE2_1[complete.cases(PoE2_1), ]
PoE2_2 <- subset(PoE2_2, select = -c(coln))

## counting each group/pop and keep one representative of each
pops <- subset(dt2, select = c(strain, pop))
pops2 <- merge(PoE2_2, pops, by= "strain")
pops2 <- data.frame(distinct(pops2) )

## keep only one strain of each group/pop
pops3 <- pops2 %>%
  group_by(Groups, pop) %>% 
  dplyr::summarise(strain=first(strain), count=n_distinct(strain) )

## adding snps data
final <- merge(pops3, dt2, by= c("strain", "pop") )
final$color <- ifelse(final$value == 1, "white", "red")
final2 <- final

final2$Subset <- factor(final2$strain, levels = unique(final2$strain))
final2$Groups <- factor(final2$Groups, levels = unique(final2$Groups))
final2$pop <- factor(final2$pop, levels = unique(final2$pop))
colors <- c(white = "#FFFFFF", red = "#000000")  

final2 <- subset(final2, color== "red")

#### removing Ts with similar patterns
final2 <- subset(final2, strain %notin% c("Br130", "PY5003-v2", "T4-2","Br48","G22","Z2-1") )
final2$Subset <- ifelse(final2$Subset == "B2", "HAP2 (B2)", final2$strain)

final2 <- subset(final2, strain %notin% c("B71", "BTGP-6e", "P29", "WB127") )
final2$Subset <- ifelse(final2$Subset == "BdBar", "HAP1 (BdBar)", final2$Subset)

# need to add U169 because all are 1 (white)
final2[nrow(final2) +1,] = c("U169","E1","47","1","4747321", "1","white","U169")


final2$Subset = factor(final2$Subset, levels = c(
  "HAP1 (BdBar)","U169","Br62","CD156","Ei8303","Ei88365","U229",
  "HAP2 (B2)", "B51",   
  "Ei9411","MG03","MG04","MG12","AR4") )

final2 <- subset(final2, sites < 4840000)

pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/FigS9/Fig_final_4", 
               ".pdf", sep = ''), 23, 4.44) 
print(
  ggplot(final2) +   #subset(df, Recipient == "ATCC64557") 
    geom_point(mapping = aes(x= as.numeric(sites), y=0, color=color, size = 2), size = 2)+  # , alpha = 0.1, show.legend = F
    scale_color_manual(name="", values=colors)+
    labs(title = " ", y = " ", x= "\nChromosome position (bp)",  colour = "Donor")+
    facet_grid(Subset ~ ., scales = "free", space = "free", switch ="y")+   #cols = vars(Chr)
    theme( plot.margin=unit(c(0, 0.8, 0, 0), units="line"), # t, r, b, l Dimensions of each margin 
           legend.title=element_text(face="bold", size=10), 
           legend.text=element_text(face="bold", size=15),
           legend.direction = "vertical",
           legend.key = element_rect(fill = "white"),
           legend.position = "right", 
           legend.justification="right",
           #legend.margin=margin(0,0,0,0),
           legend.box.margin=margin(0,0,0,0.1, "cm"),
           legend.spacing.y = unit(0.01, "cm"), 
           axis.title.x = element_text(face="bold", size = 17),  # , margin = margin(1, 0, 40, 0), vjust = -3
           axis.text=element_text(face="bold", angle=0,  size = 12, colour = "black"), #element_text(face="bold", angle=0,  size = 8, colour = "black"), 
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),    #element_line(size=1),
           axis.ticks.length.x = unit(0.2,"cm"),
           axis.ticks.x = element_line(size = 0.4),
           strip.text.x = element_text(size=12, face = "bold", colour = "black", margin = margin(0.15,0,0.15,0, "cm")),
           strip.text.y.left = element_text(size=13, face = "bold", angle = 0, colour = "black", margin = margin(0.15,0.15,0.15,0.15, "cm")),
           strip.background=element_rect(fill = "grey", colour = "black", size = 0.35),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(), 
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(), #element_line(colour="gray", size=0.5),
           panel.background = element_blank(), #element_rect(fill = "white", colour = "black", size = 0.2),
           panel.border = element_rect(fill = NA, colour = "gray", size = 0.35),#element_rect(colour ="white", fill=NA, size=5) 
           panel.spacing.x=unit(0.2, "lines"),
           panel.spacing.y=unit(0.1, "lines"),
           axis.line.x = element_line(colour = "black", size = 0.1, linetype = "solid")
    )
  + guides(colour = guide_legend(override.aes = list(size=2.5)))
  
)
dev.off()

######################## PoSt
'%notin%' <- Negate('%in%')

dtst <- read.csv("~/Google Drive/1.NATURE_GENETICS/FigS9/PoSt/2PoSt_Chr7_1790_2000.csv", header=TRUE)
dtst <- subset(dtst, select = -c(col))
dtst1 <- subset(dtst, pop %in% c("T", "Lee", "P", "O", "S", "P2") ) #   
#dtst1$pop <- ifelse(dtst1$pop=="Lee"| dtst1$pop=="P"|dtst1$pop=="O"|dtst1$pop=="S"|dtst1$pop=="P2", "X", dtst1$pop)

#choosing the first 200 columns
dtst1 <- dtst1 %>% select(1:220)
dt2 <- gather(dtst1, sites, value, X171523:X177915)
dt2$sites <- gsub("X", "", as.character(dt2$sites))
###### choosing only St s
dtst <- read.csv("~/Google Drive/1.NATURE_GENETICS/FigS9/PoSt/2PoSt_Chr7_1790_2000.csv", header=TRUE)
dtst <- subset(dtst, select = -c(col))
dtst1 <- subset(dtst, pop %in% c("St") ) #   
#choosing the first 200 columns
dtst1 <- dtst1 %>% select(1:220)
dt2_st <- gather(dtst1, sites, value, X171523:X177915)
dt2_st$sites <- gsub("X", "", as.character(dt2_st$sites))

################################# choosing only E s
dtst <- read.csv("~/Google Drive/1.NATURE_GENETICS/FigS9/PoSt/2PoSt_Chr7_1790_2000.csv", header=TRUE)
dtst <- subset(dtst, select = -c(col))
dtst1 <- subset(dtst, pop %in% c("T", "E1", "E2", "E3") ) #   
dtst1$pop <- ifelse(dtst1$pop=="E1"| dtst1$pop=="E2"|dtst1$pop=="E3", "E", dtst1$pop)
#choosing the first 200 columns
dtst1 <- dtst1 %>% select(1:220)
dt2_ET <- gather(dtst1, sites, value, X171523:X177915)
dt2_ET$sites <- gsub("X", "", as.character(dt2_ET$sites))

############################################################
### reading grouping data
PoSt <- read_excel("~/Google Drive/1.NATURE_GENETICS/FigS9/PoSt/New_hapl_grps.xlsx")
PoSt_1 <- gather(PoSt, coln, strain, V1:V58)
PoSt_2 <- PoSt_1[complete.cases(PoSt_1), ]
PoSt_2 <- subset(PoSt_2, select = -c(coln))

## counting each group/pop and keep one representative of each
pops <- subset(dt2, select = c(strain, pop))
pops2 <- merge(PoSt_2, pops, by= "strain")
pops2 <- data.frame(distinct(pops2) )

# group_by(Groups, pop) %>% dplyr::summarise(strain, count= n()) 
## keep only one strain of each group/pop
pops3 <- pops2 %>%
  group_by(Groups, pop) %>% 
  dplyr::summarise(strain=first(strain), count=n_distinct(strain) )
## adding snps data
final_X <- merge(pops3, dt2, by= c("strain", "pop") )
final_X <- subset(final_X, select = -c(count))
final_X$pop <- ifelse(final_X$pop=="Lee"| final_X$pop=="P"|final_X$pop=="O"|final_X$pop=="S"|final_X$pop=="P2", "X", final_X$pop)
################################ making st
pops <- subset(dt2_st, select = c(strain, pop))
pops2 <- merge(PoSt_2, pops, by= "strain")
pops2 <- data.frame(distinct(pops2) )
final_st <- merge(pops2, dt2_st, by= c("strain", "pop") )
############################### Es
pops <- subset(dt2_ET, select = c(strain, pop))
pops2 <- merge(PoSt_2, pops, by= "strain")
pops2 <- data.frame(distinct(pops2) )

pops3 <- pops2 %>%
  group_by(Groups, pop) %>% 
  dplyr::summarise(strain=first(strain), count=n_distinct(strain) )
## adding snps data
final_ET <- merge(pops3, dt2_ET, by= c("strain", "pop") )
final_ET <- subset(final_ET, select = -c(count))


#############################
Final <- rbind(final_st,final_X,final_ET)

Final$Subset <- Final$strain
#### removing STs with similar patterns
Final <- subset(Final, strain %notin% c("SSFL02-1","SSFL14-3","SSTX16-11","STAG-MS","U217","U233") )
Final$Subset <- ifelse(Final$Subset == "Br35", "PoST/U3 (Br35)", Final$strain)
Final$Subset <- ifelse(Final$Subset == "Pg1054", "PoST (Pg1054)", Final$Subset)
Final$Subset <- ifelse(Final$Subset == "UbJA112", "PoST (UbJA112)", Final$Subset)

############################## T and X s
Final <- subset(Final, strain %notin% c("Arcadia2","Ken53-33","87-120","Bd8401","BP1","YNPM4-1-1","BTTrp-5","GFSI1-7-2","IA1","IB33",
                                        "Lh8401","Lh88405-2","MG05","MG08","ML33","NNPM1-2-1", "Sv9610", "U232",
                                        "Ei9411", "AR4") )
Final$Subset <- ifelse(Final$Subset == "Br130", "PoX (Br130-T)", Final$Subset)
Final$Subset <- ifelse(Final$Subset == "Br108-1", "PoSt/U3 (Br108-1-T)", Final$Subset)
#Final$Subset <- ifelse(Final$Subset == "B71", "PoE (B71-T)", Final$Subset)
Final$Subset <- ifelse(Final$Subset == "B2", "PoSt (B2-T)", Final$Subset)
Final$Subset <- ifelse(Final$Subset == "Sv9623", "PoX (Sv9623)", Final$Subset)

Final <- subset(Final, value== "0")
Final$color <- ifelse(Final$pop == "St" , "St", "other")
Final$color <- ifelse(Final$pop == "X" , "X", Final$color)
Final$color <- ifelse(Final$Subset == "PoX (Br130-T)" , "X", Final$color)
Final$color <- ifelse(Final$Subset == "PoSt/U3 (Br108-1-T)" , "St", Final$color)
Final$color <- ifelse(Final$Subset == "PoSt (B2-T)" , "St", Final$color)
Final$color <- ifelse(Final$pop == "E" , "E", Final$color)
Final$color <- ifelse(Final$value == "1" , "one", Final$color)

Final[nrow(Final) +1,] = c("B71","T","28","176747", "1","PoE (B71-T)", "one")
Final[nrow(Final) +1,] = c("B51","T","28","172000", "1","PoE (B51-E1)", "one")

colors <- c(St= "#000000",  X= "#000000", E = "#000000", one = "#FFFFFF") 
Final$Subset = factor(Final$Subset, levels = c("PoSt/U3 (Br108-1-T)",
                                                "PoST/U3 (Br35)",
                                               "PoSt (B2-T)",
                                                "PoST (UbJA112)",
                                               "PoST (Pg1054)", ####"PoT (B2)","PoST (SSFL02-1)",  "PoST (SSFL14-3)", "PoST (STAG-MS)","PoST (U217)","PoST (U233)","PoST (SSTX16-11) -3",
                                                "PoX (Br130-T)", "PoX (Sv9623)",  
                                                "PoE (B71-T)", "PoE (B51-E1)"))#

pdf(file=paste("~/Google Drive/1.NATURE_GENETICS/FigS9/PoSt/Fig_final_4", 
                              ".pdf", sep = ''), 23, 4) 
print(
  ggplot(Final) +   #subset(df, Recipient == "ATCC64557") 
    geom_point(mapping = aes(x= as.numeric(sites), y=0, color=color, size = 2), size = 2)+  # , alpha = 0.1, show.legend = F
    scale_color_manual(name="", values=colors)+
    labs(title = " ", y = " ", x= "\nChromosome position (bp)",  colour = "Donor")+
    facet_grid(Subset ~ ., scales = "free", space = "free", switch ="y")+   #cols = vars(Chr)
    theme( plot.margin=unit(c(0, 0.8, 0, 0), units="line"), # t, r, b, l Dimensions of each margin 
           legend.title=element_text(face="bold", size=10), 
           legend.text=element_text(face="bold", size=15),
           legend.direction = "vertical",
           legend.key = element_rect(fill = "white"),
           legend.position = "right", 
           legend.justification="right",
           #legend.margin=margin(0,0,0,0),
           legend.box.margin=margin(0,0,0,0.1, "cm"),
           legend.spacing.y = unit(0.01, "cm"), 
           axis.title.x = element_text(face="bold", size = 17),  # , margin = margin(1, 0, 40, 0), vjust = -3
           axis.text=element_text(face="bold", angle=0,  size = 12, colour = "black"), #element_text(face="bold", angle=0,  size = 8, colour = "black"), 
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),    #element_line(size=1),
           axis.ticks.length.x = unit(0.2,"cm"),
           axis.ticks.x = element_line(size = 0.4),
           strip.text.x = element_text(size=12, face = "bold", colour = "black", margin = margin(0.15,0,0.15,0, "cm")),
           strip.text.y.left = element_text(size=13, face = "bold", angle = 0, colour = "black", margin = margin(0.15,0.15,0.15,0.15, "cm")),
           strip.background=element_rect(fill = "grey", colour = "black", size = 0.35),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(), 
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(), #element_line(colour="gray", size=0.5),
           panel.background = element_blank(), #element_rect(fill = "white", colour = "black", size = 0.2),
           panel.border = element_rect(fill = NA, colour = "gray", size = 0.35),#element_rect(colour ="white", fill=NA, size=5) 
           panel.spacing.x=unit(0.2, "lines"),
           panel.spacing.y=unit(0.1, "lines"),
           axis.line.x = element_line(colour = "black", size = 0.1, linetype = "solid")
    )
  + guides(colour = guide_legend(override.aes = list(size=2.5)))
  
)
dev.off()
