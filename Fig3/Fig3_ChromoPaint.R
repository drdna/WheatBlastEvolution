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

#strains <- read_excel("~/Google Drive/1.SCIENCE_WHEAT_BLAST/genotypesFile_Ed copy.xlsx", sheet = "newSet")
strains <- read_excel("~/Desktop/AUG2022/REVISION2/HaplotypesTable.xlsx", sheet = "Sheet1")

#strains[is.na(strains)] <- ""
strains <- select(strains, c(2,3,15))
#strains$geneo <- paste(strains$Host, strains$genotype, sep = " ")




#strains <- subset(strains, (Geno %in% c( "PoL1-1","PoL1-1","PoL1-1","PoL1-1","PoL1-1","PoL1-1",
#                                         "PoL1-1","PoT1","PoL1-2", "PoL1-6","PoL1-12","PoL1-13", "PoT1", "PoT4", "PoT32"))) #,
#"PoL1-9","PoL1-11","PoL1-12", "PoL1-13", "PoT2","PoT3","PoT4","PoT5","PoT6",
#"PoT7","PoT8","PoT9","PoT10","PoT11", "PoT12","PoT13","PoT14","PoT15","PoT16",
#"PoT17","PoT18","PoT19","PoT20","PoT21", "PoT23","PoT24","PoT25","PoT26","PoT27",
#"PoT29","PoT30","PoT31","PoT32", "PoT34")))



strains <- subset(strains, (Geno %in% c( "PoL1-1","PoT1","PoL1-2","PoL1-3","PoL1-4","PoL1-5","PoL1-6","PoL1-7","PoL1-8",
                                         "PoL1-9","PoL1-11","PoL1-12", "PoL1-13", "PoT2","PoT3","PoT4","PoT5","PoT6",
                                         "PoT7","PoT8","PoT9","PoT10","PoT11", "PoT12","PoT13","PoT14","PoT15","PoT16",
                                         "PoT17","PoT18","PoT19","PoT20","PoT21", "PoT23","PoT24","PoT25","PoT26","PoT27",
                                         "PoT29","PoT30","PoT31","PoT32", "PoT34")))

file.dir <- "~/CPtest/"
setwd(file.dir)
lengths<- list(6442611, 7903502, 8206720, 5402944, 4443078, 6091895, 4042924)

#Chromosomes <- list("chr2")
Chromosomes <- list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7")

count <- 0
tabless <- list()
for(strain in strains$Recipient) {
  tables <- list()
  for (chr in Chromosomes ){
    copyprobsfile = paste0(strain,".",chr,".V2.complete_out.forR.copyprobsperlocus.out")
        count <- as.numeric(count + 1) 
        system2("echo", shQuote(lengths[[count]]))
        system2("echo", shQuote(count))
  
        # read in copyprobsfile
        dt <- read.table(paste(file.dir, copyprobsfile, sep=""), header=TRUE)
         # adding E1 and E2 to E
        dt <- dt %>% mutate(E = E1+E2+E3) 
        dt <- dt %>% mutate(X = Lee+P+O+S+P2)
        dt <- subset(dt, select= c("pos",  "E","Lu", "St", "U1","X"))
        ################################### choosing the population with the highest copyprobs 2times of the next highest one
        process_row = function(line) {
            firstMax = max(line)
            secMax = sort(line, TRUE)[2]
            # secMax2 = max(line[line != firstMax])
              #print(firstMax)
                #print(secMax)
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
      #identify corresponding widths file
      widthsfile <- (gsub("\\.forR.copyprobsperlocus.out$","\\.width.copyprobsperlocus.out", copyprobsfile))
      system2("echo", shQuote(widthsfile))
  
      #Read in bar width vector
      width <- t(read.table(widthsfile, header=FALSE, row.names=1))
      width <- data.frame(width)
      ###condition to play with width size
      avr = mean(width$width)+20000
      width$width2 <- ifelse(width$width > avr , avr, width$width)
      #Grab values in "center" column of widthsfile as column headers for transposed matrix
      names <-  width$center # (width[,c('center')])
  
      #transpose data matrix after removing column 1 (positions)
      dtt <- as.data.frame(t(dt[,-1]))
  
       # add column names to transposed matrix
        colnames(dtt) <- names
  
      #bind columns needed for plotting
      test <- cbind(dtt, donor = rownames(dtt))
  
      #melt data for plotting
      dm = reshape2::melt(test, id.vars = c('donor'))
  
      # add width values to melted matrix
      #dm$width <- rep((width$width2), each=ncol(dt)-1)
  
      # grab every 2000th column label for x-axis labels
      labels <- as.data.frame(names)
      #Required: convert SNP positions to numeric values
      dm$variable=as.numeric(levels(dm$variable))[dm$variable]
      #reverse the matrix
      reversed_dm <- dplyr::arrange(dm, -row_number())
  
      ## getting chr name

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
#df[df==0] <- NA
#adding genotype names based on strains table
df2 <- as.data.frame(merge(df, strains, by = "Recipient"))

df2$StrainDate <- paste(df2$Geno, " (", df2$date, ")", sep = "")

#df2[str_sort(df2$StrainDate, numeric = TRUE), ]

df2 <- subset(df2, select = -c(Recipient, date))

#Changing Chr to Chromosome
df2$Chr <- gsub("Chr", "Chromosome ", df2$Chr)


df2 <- as.data.frame(df2)

#df2 <- read.csv("~/Google Drive/1.NATURE_GENETICS/Fig4/Final Figure/Others/DF_for_Fig4_Final_Final.csv")
df2 <- merge(df2, strains, by = "Geno")

df2$Subset = factor(df2$StrainDate, #levels = unique(str_sort(df2$StrainDate, numeric = TRUE)))
levels = c("PoL1-1 (1980)",  "PoT1 (1985)", "PoL1-2 (1988)",  "PoL1-3 (1990)",
           "PoL1-4 (2002)",  "PoL1-5 (2002)",   "PoL1-6 (2005)",  "PoL1-7 (2008)",
           "PoL1-8 (2010)",  "PoL1-9 (2012)",  "PoL1-11 (2014)", "PoL1-12 (2014)",
           "PoL1-13 (2017)", "PoT2 (1986)",    "PoT3 (1987)",    "PoT4 (1987)",
           "PoT5 (1988)",    "PoT6 (1988)",    "PoT7 (1989)",    "PoT8 (1989)",
           "PoT9 (1989)",    "PoT10 (1990)",   "PoT11 (1990)",   "PoT12 (1990)",
           "PoT13 (1991)",   "PoT14 (1992)",   "PoT15 (1992)",   "PoT16 (1992)",
           "PoT17 (1992)",   "PoT18 (1992)",   "PoT19 (1992)",   "PoT20 (1992)",
           "PoT21 (1992)",   "PoT23 (1992)",   "PoT24 (1992)",   "PoT25 (1992)",
           "PoT26 (2005)",   "PoT27 (2005)",   "PoT29 (2007)",   "PoT30 (2009)",
           "PoT31 (2012)",   "PoT32 (2012)",   "PoT34 (2012)"))

#oldColors <- c( 'E' = "#A3A500",Lu= "#CC3333", St= "#F564E3", U1 = "#F8766D", X= "#619CFF")  #, E3="#717200"

colors <- c( 'E' = "#009E73",Lu= "#999999", St= "#AA4499", U1 = "#D55E00", X= "#0072B2") 

ppi <- 300
#tiff(file=paste("~/Google Drive/1.NATURE_GENETICS/Fig4/Final Figure/Fig4_Final_Final_withStrain",".tiff", sep = ''), 
#     width = 14*ppi, height= 7.8*ppi, res = ppi) 
pdf(file=paste("~/Fig3-1",".pdf", sep = ''), width = 20, height= 8.5) 

# subsample dataset every fifth line
#df2$variable <- df2$variable/1000000

#df2$variable <- round(df2$variable, 1)

#df2 <- df2[!duplicated(c(df2$variable)), ]

print(
  ggplot(df2) +   #subset(df, Recipient == "ATCC64557") 
    geom_point(mapping = aes(x=variable, y=0, color=donor, size = 3.5), size = 3.5)+  # , alpha = 0.1, show.legend = F
    scale_color_manual(name="Donor", values=colors)+
    scale_x_continuous(minor_breaks = seq(0, 8, 0.1) )+
    coord_cartesian(ylim = c(-0.5,0.5)) + 
    labs(title = " ", y = " ", x= "\nChromosome position (Mb)",  colour = "Donor")+
    facet_grid(Subset ~ Chr, scales = "free", space = "free", switch ="y")+   #cols = vars(Chr)
    
    theme( 
           plot.margin=unit(c(0.1, 0.8, 0.1, 0.1), units="line"), # t, r, b, l Dimensions of each margin 
           legend.title=element_text(face="bold", size=10), 
           legend.text=element_text(face="bold", size=9),
           legend.direction = "vertical",
           legend.key = element_rect(fill = "white"),
           legend.position = "right", 
           legend.justification="right",
           #legend.margin=margin(0,0,0,0),
           legend.box.margin=margin(0,0,0,0.1, "cm"),
           legend.spacing.y = unit(0.01, "cm"), 
           axis.title.x = element_text(face="bold", size = 13),  # , margin = margin(1, 0, 40, 0), vjust = -3
           axis.text=element_text(face="bold", angle=0,  size = 8, colour = "black"), 
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),    #element_line(size=1),
           axis.ticks.length.x = unit(0.15,"cm"),
           axis.ticks.x = element_line(size = 0.3),
           strip.text.x = element_text(size=8, face = "bold", colour = "black", margin = margin(0.15,0,0.15,0, "cm")),
           strip.text.y.left = element_text(size=9, face = "bold", angle = 0, colour = "black", margin = margin(0.15,0.15,0.15,0.15, "cm")),
           strip.background=element_rect(fill = "grey", colour = "black", size = 0.35),
           strip.placement = "outside",
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(), 
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(), #element_line(colour="gray", size=0.5),
           panel.background = element_blank(), #element_rect(fill = "white", colour = "black", size = 0.2),
           panel.border = element_blank(),#element_rect(colour ="white", fill=NA, size=5) 
           panel.spacing.x=unit(0.2, "lines"),
           panel.spacing.y=unit(0.1, "lines"),
           axis.line.x = element_line(colour = "black", size = 0.1, linetype = "solid")
    )
  + guides(colour = guide_legend(override.aes = list(size=5)))
)

dev.off()




