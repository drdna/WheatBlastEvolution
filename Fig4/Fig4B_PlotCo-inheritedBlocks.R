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

# Read table to select a single reference strain for each haplotype (first column)
strains <- read_excel("/Users/mfarman/Desktop/AUG2022/REVISION2/HaplotypesTable.xlsx", sheet = "Sheet1")

#only grab columns 2, 3 & 15
strains <- select(strains, c(2,3,15))

# list strains to include in analysis (poor quality genomes dropped)
strains <- subset(strains, (Geno %in% c( "PoL1-1","PoT1","PoL1-2","PoL1-3","PoL1-4","PoL1-5","PoL1-6",
                                         "PoL1-7","PoL1-8","PoL1-9",
                                         "PoL1-11","PoL1-12","PoL1-13",  
                                           "PoT2","PoT3","PoT4","PoT5","PoT6","PoT7","PoT8","PoT9","PoT10","PoT11",
                                           "PoT12","PoT13","PoT14","PoT15","PoT16",
                                         "PoT17","PoT18","PoT19","PoT20","PoT21",
                                         "PoT23", "PoT24","PoT25","PoT26","PoT27",
                                         "PoT29","PoT30","PoT31","PoT32",                "PoT34") ))

# path to ChromoPainter (CP) outfiles
file.dir <- "~/CPtest/"
setwd(file.dir)

# chromosome lengths
lengths<- list(6442611, 7903502, 8206720, 5402944, 4443078, 6091895, 4042924)

# chromosome names in CP filenames
Chromosomes <- list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7")

# Chromosome lengths stored in dataframe (to allow different x-axes for facets)
ChrLengths = data_frame(Subset = NA, length = rep(c(6.442611, 7.903502, 8.206720, 5.402944, 4.443078, 6.091895, 4.042924), each = 4),
                        Chr = rep(c("Chromosome 1", "Chromosome 2", "Chromosome 3", "Chromosome 4", "Chromosome 5", "Chromosome 6", "Chromosome 7"), each = 4))

# Add subset info to dataframe for plotting different groups of haplotypes
ChrLengths$Subset <- rep(c("PoL", "Lolium", "PoT", "Wheat"), 7)


# Use loop to read in copyprobs files

count <- 0

tabless <- list()
for(strain in strains$Recipient) {
  tables <- list()
  for (chr in Chromosomes ){
    #strain="CHRF"
    #chr="chr2"
    copyprobsfile = paste0(strain,".",chr,".V2.complete_out.forR.copyprobsperlocus.out")
        #count the loops
        count <- as.numeric(count + 1) 
        system2("echo", shQuote(lengths[[count]]))
        system2("echo", shQuote(count))
  
        # read copyprobsfile
        dt <- read.table(paste(file.dir, copyprobsfile, sep=""), header=TRUE)
         # adding E1 and E2 to E
        dt <- dt %>% mutate(E = E1+ E2+E3) 
        dt <- dt %>% mutate(X = Lee+P+O+S+P2)
        dt <- subset(dt, select= c("pos",  "E","Lu", "St", "U1","X"))
        
        ################################### choosing the population with the highest copyprobs ≥2 times the next highest one
        process_row = function(line) {
            firstMax = max(line)
            secMax = sort(line, TRUE)[2]
                if (firstMax > 2 * secMax) {
                  line[line != firstMax] = 0
                        #print(line)
                          }
                    else{
                        line[] = 0
                                            }
                          return(line)
                    }
        ###########
        x <- dt[,-1]
        new_dt = apply(x, MARGIN = 1, FUN = process_row)
        dt2 = as.data.frame(t(new_dt))
        dt <- as.data.frame(cbind(dt[1], dt2))
        ###########################################################################################################################################
  
      #transpose data matrix after removing column 1 (positions)
      positions <- as.numeric(dt[,1])
      dtt <- as.data.frame(t(dt[,-1]))

      #bind columns needed for plotting
      test <- cbind(dtt, donor = rownames(dtt))
  
      #melt data for plotting
      dm = reshape2::melt(test, id.vars = c('donor'), variable.name = 'index')
      
      dm$Pos <- as.numeric(positions[dm$index])
      dm$Pos=as.numeric(dm$Pos)
      
      #reverse the matrix (CP gives results in reverse order)
      reversed_dm <- dplyr::arrange(dm, -row_number())
  
      # getting chr name

      name = paste0(str_sub(chr, 2, -1))
      name = paste("C", name, sep = "")

      reversed_dm$Chr <- name
  
      reversed_dm$Recipient <- strain
  
      tables[[name]]= reversed_dm
  }
  
  dddd=data.table::rbindlist(tables, fill = T)
  
  table_name = paste(strain, "Table", sep = "_")
  
  tabless[[strain]] = dddd
  count <- 0 
}  

d_circle=data.table::rbindlist(tabless, fill = T)
  
########### make circles for 7 chrs for each strains
df <- subset(d_circle, value != 0) 

# uncomment next line to generate an output file from dataframe
#write.csv(df,"~/Google Drive/1.NATURE_GENETICS/Fig4B/df.csv", row.names = FALSE)

# add haplotype IDs based on strains table
df2 <- merge(df, strains, by = "Recipient")
df2 <- as.data.frame(df2)
df2 <- subset(df2, select = -c(Recipient, date))

# Changing Chr to Chromosome
df2$Chr <- gsub("Chr", "Chromosome ", df2$Chr)

# specify PoT isolate subset
df3 <- subset(df2, (Geno %in% c("PoT1","PoT2","PoT3","PoT4","PoT5","PoT6","PoT7","PoT8","PoT9","PoT10","PoT11",
                                "PoT12","PoT13","PoT14","PoT15","PoT16",
                                "PoT17","PoT18","PoT19","PoT20","PoT21",
                                "PoT23", "PoT24","PoT25","PoT26","PoT27",
                                "PoT29","PoT30","PoT31","PoT32",                "PoT34") ))


## identify positions where all PoT members have same donor
df3$sameDonor <- df3$donor
df4 <- df3 %>%
  group_by(Pos, Chr) %>%
  mutate(sameDonor = replace(sameDonor, n_distinct(sameDonor)==1, 'True') )

df4 <- subset(df4, sameDonor == 'True')

df5 <- df4 %>%
  group_by(Pos, Chr) %>% summarise(donor, Geno, count= n(), .groups = 'drop') 
df5 <- subset(df5, count ==max(unique(df5$count)))

df_T <- subset(df5, (Geno %in% c( "PoT1") ))
df_T$Geno <- "PoT"



# specify PoL isolate subset 
df_L <- subset(df2, (Geno %in% c("PoL1-1","PoL1-2","PoL1-3","PoL1-4","PoL1-5","PoL1-6",
                                 "PoL1-7","PoL1-8","PoL1-9",
                                 "PoL1-11","PoL1-12","PoL1-13"  )))

## identify positions where all PoL members have same donor
df_L$sameDonor <- df_L$donor
df_L1 <- df_L %>%
  group_by(Pos, Chr) %>%
  mutate(sameDonor = replace(sameDonor, n_distinct(sameDonor)==1, 'True') )

df_L1 <- subset(df_L1, sameDonor == 'True')

df_L2 <- df_L1 %>%
  group_by(Pos, Chr) %>% summarise(donor, Geno, count= n(), .groups = 'drop') 
df_L2 <- subset(df_L2, count ==max(unique(df_L2$count)))
df_L3 <- subset(df_L2, (Geno %in% c( "PoL1-1") ))
df_L3$Geno <- "PoL"

# specify all haplotypes that were found on wheat
df_Lw <- subset(df2, (Geno %in% c("PoL1-5", "PoL1-6", "PoL1-7", 
                                  "PoT2","PoT3","PoT4","PoT5","PoT6","PoT7","PoT8","PoT9","PoT10","PoT11",
                                  "PoT12","PoT13","PoT14","PoT15","PoT16",
                                  "PoT17","PoT18","PoT19","PoT20","PoT21",
                                  "PoT23", "PoT24","PoT25","PoT26","PoT27",
                                  "PoT29","PoT30","PoT31","PoT32",                "PoT34")))



## identify positions where all isolates found on wheat members have same donor
df_Lw$sameDonor <- df_Lw$donor
df_Lw1 <- df_Lw %>%
  group_by(Pos, Chr) %>%
  mutate(sameDonor = replace(sameDonor, n_distinct(sameDonor)==1, 'True') )

df_Lw1 <- subset(df_Lw1, sameDonor == 'True')

df_Lw2 <- df_Lw1 %>%
  group_by(Pos, Chr) %>% summarise(donor, Geno, count= n(), .groups = 'drop') 
df_Lw2 <- subset(df_Lw2, count ==max(unique(df_Lw2$count)))
df_Lw3 <- subset(df_Lw2, (Geno %in% c( "PoL1-5") ))
df_Lw3$Geno <- "Wheat"
  
# specify haplotypes that were found on Lolium
df_Lolium <- subset(df2, (Geno %in% c ("PoL1-1","PoL1-3","PoL1-4","PoL1-5","PoL1-6",
                                         "PoL1-8",      "PoL1-13"  
                                       ) ))

## identify positions where all isolates found on Lolium have same donor
df_Lolium$sameDonor <- df_Lolium$donor
df_Lolium1 <- df_Lolium %>%
  group_by(Pos, Chr) %>%
  mutate(sameDonor = replace(sameDonor, n_distinct(sameDonor)==1, 'True') )

df_Lolium1 <- subset(df_Lolium1, sameDonor == 'True')

df_Lolium2 <- df_Lolium1 %>%
  group_by(Pos, Chr) %>% summarise(donor, Geno, count= n(), .groups = 'drop') 
df_Lolium2 <- subset(df_Lolium2, count ==max(unique(df_Lolium2$count)))
df_Lolium3 <- subset(df_Lolium2, (Geno %in% c( "PoL1-1") ))
df_Lolium3$Geno <- "Lolium"


# creat a final datafarme for plotting

df_final <- rbind(df_T, df_L3, df_Lw3,df_Lolium3)
df_final$Subset = factor(df_final$Geno, levels = c( "PoL", "Lolium","PoT","Wheat") )

# old color scheme
#colors <- c( 'E' = "#A3A500",Lu= "#CC3333", St= "#F564E3", U1 = "#F8766D", X= "#619CFF")  #, E3="#717200"

# current color scheme
colors <- c( 'E' = "#009E73",Lu= "#999999", St= "#AA4499", U1 = "#D55E00", X= "#0072B2")

ppi <- 300

#uncomment to generate output file
#tiff(file=paste("~/Google Drive/Farman_Lab/1.NATURE_GENETICS/1.SUBMISSION/1.Revison_2/Final_Fig4",".tiff", sep = ''), 
#     width = 13*ppi, height= 2.4*ppi, res = ppi) #height= 7.8*ppi,
#pdf(file=paste("~/Final_Fig4B",".pdf", sep = ''), 20, 2.4)


p <- ggplot(df_final) +   # subset(df, Recipient == "ATCC64557") # uncomment to plot only certain strains
    geom_point(mapping = aes(x=Pos, y=0, color=donor, size = 3.5), size = 3.5) + 
    scale_color_manual(name="Donor", values=colors) +
    scale_x_continuous(breaks = seq(0, 8, 2)) + 
    coord_cartesian(ylim = c(-0.5,0.5)) + 
    labs(title = " ", y = " ", x= "\nChromosome position (Mb)",  colour = "Donor")+
    facet_grid(Subset ~ Chr, scales = "free", space = "free", switch ="y")+ 
    theme( plot.margin=unit(c(0, 0.8, 0, 0), units="line"), # t, r, b, l Dimensions of each margin 
           legend.title=element_text(face="bold", size=10), 
           legend.text=element_text(face="bold", size=9),
           legend.direction = "vertical",
           legend.key = element_rect(fill = "white"),
           legend.position = "right", 
           legend.justification="right",
           legend.box.margin=margin(0,0,0,0.1, "cm"),
           legend.spacing.y = unit(0.01, "cm"), 
           axis.title.x = element_text(face="bold", size = 13),
           axis.text=element_text(face="bold", angle=0,  size = 8, colour = "black"), 
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.ticks.length.x = unit(0.15,"cm"),
           axis.ticks.x = element_line(size = 0.3),
           strip.text.x = element_text(size=8, face = "bold", colour = "black", margin = margin(0.15,0,0.15,0, "cm")),
           strip.text.y.left = element_text(size=9, face = "bold", angle = 0, colour = "black", margin = margin(0.15,0.15,0.15,0.15, "cm")),
           strip.background=element_rect(fill = "grey", colour = "black", size = 0.35),
           strip.placement = "outside",
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(), 
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(colour ="grey", fill=NA, size=1),
           panel.spacing.x=unit(0.2, "lines"),
           panel.spacing.y=unit(0.1, "lines"),
           axis.line.x = element_line(colour = "black", size = 0.1, linetype = "solid")
    ) + guides(colour = guide_legend(override.aes = list(size=5))) + expand_limits(x = 0)
p
p1 <- p + geom_blank(data = ChrLengths, aes(x=length, y=0)) +
  facet_grid(Subset ~ Chr, scales= "free", space = "free", switch = "y")

p1
#dev.off()


