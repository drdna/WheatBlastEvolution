
#Ploting ChroPainter's copyprobsperlocus data for having both circle (dominant population) and barplots(all population)
#Author: Mostafa Rahnama 

#open libraries
#Rcode for plotting chromopaintings along continuous x-axis
library(reshape2)
library(scales)
library(ggplot2)
library(grid)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)


Chromosomes <- list("Chromosome 1", "Chromosome 2", "Chromosome 3", "Chromosome 4", "Chromosome 5", "Chromosome 6", "Chromosome 7")
lengths<- list(6442611, 7903502, 8206720, 5402944, 4443078, 6091895, 4042924)

############ Reading data for bar plots 
file.dir <- "~/Out_3/"
setwd("~/Out_3/")

count <- 0
tables <- list()
for(copyprobsfile in dir(file.dir, pattern="^B71.*\\V2.complete_out.forR.copyprobsperlocus.out$")) {    #pattern="^B71.*\\_Paint_strainsID10.copyprobsperlocus.out.forR$"
  #count the loops                         
  count <- as.numeric(count + 1)
  system2("echo", shQuote(lengths[[count]]))
  system2("echo", shQuote(count))
  # read in copyprobsfile
  #copyprobsfile <- "/Users/mostafa/Google Drive/1.NATURE_GENETICS/CPApr042921/Out_3/B71.chr1.V2.complete_out.forR.copyprobsperlocus.out"
  dt <- read.table(copyprobsfile, header=TRUE)
  # adding E1 and E2 to E
  dt <- dt %>% mutate(E = E1+ E2) 
  ## remove E1 and E2
  dt <- subset(dt, select = -c(E1, E2))
  colnames(dt)[15] <- "E1+E2"
  #identify corresponding widths file
  widthsfile <- (gsub("\\.forR.copyprobsperlocus.out$","\\.width.copyprobsperlocus.out", copyprobsfile))
  system2("echo", shQuote(widthsfile))
  #Read in bar width vector
  width <- t(read.table(widthsfile, header=FALSE, row.names=1))
  width <- data.frame(width)
  ###condition to play with width size
  avr = mean(width$width)+20000
  #if width greater than avr replace witdth with avr
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
  dm$width <- rep((width$width2), each=ncol(dt)-1)
  # grab every 2000th column label for x-axis labels
  labels <- as.data.frame(names)
  #Required: convert SNP positions to numeric values
  dm$variable=as.numeric(levels(dm$variable))[dm$variable]
  #reverse the matrix
  reversed_dm <- dplyr::arrange(dm, -row_number())
  name = paste0(substr(copyprobsfile, start = 6, stop = 8))
  name = paste("C", name, sep = "")
  reversed_dm$Chr <- name
  
  tables[[name]]= reversed_dm
  
}  
d=data.table::rbindlist(tables, fill = T)


############ Reading data for circle plots 

count <- 0
tables <- list()
for(copyprobsfile in dir(file.dir, pattern="^B71.*\\V2.complete_out.forR.copyprobsperlocus.out$")) {
  
  #count the loops
  count <- as.numeric(count + 1)
  system2("echo", shQuote(lengths[[count]]))
  system2("echo", shQuote(count))
  # read in copyprobsfile
  dt <- read.table(copyprobsfile, header=TRUE)
  # adding E1 and E2 to E
  dt <- dt %>% mutate(E = E1+ E2) 
  ## remove E1 and E2
  dt <- subset(dt, select = -c(E1, E2))
  colnames(dt)[15] <- "E1+E2"
  ################################### chosing the population with the highest copyprobs 2times of the next highest one
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
  #################################################################
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
  dm$width <- rep((width$width2), each=ncol(dt)-1)
  # grab every 2000th column label for x-axis labels
  labels <- as.data.frame(names)
  #Required: convert SNP positions to numeric values
  dm$variable=as.numeric(levels(dm$variable))[dm$variable]
  #reverse the matrix
  reversed_dm <- dplyr::arrange(dm, -row_number())
  
  name = paste0(substr(copyprobsfile, start = 6, stop = 8))
  name = paste("C", name, sep = "")
  reversed_dm$Chr <- name
  tables[[name]]= reversed_dm
}  
d_circle=data.table::rbindlist(tables, fill = T)

#########################################################################################
colors <- c(Er="#CC6666",  Lu= "#999999", 
            Ec="#669999",  E3="#717200",
            U2 = "#D2B48C", Lee="#A3A500", O="#00CCCC",  
            S= "#0072B2", 'E1+E2' = "#009E73", U1 = "#D55E00", P = "#00CC66", St = "#AA4499",
            C1="#FFA500", P2="#9370DB" ) #<<< new , missed=L2="#666600",P1="#9370DB", C2= "#FFA500",

colnames(d) <- c("donor", "variable", "valueD", "width", "Chr")
Merge_DT <- d_circle %>% right_join(d, by=c("donor","variable", "width", "Chr"))
Merge_DT <- Merge_DT %>% replace(is.na(.), 0)

pdf(file="~/Fig2A.pdf", 48, 40)
print(
  ggplot(Merge_DT) + 
    geom_bar(mapping = aes(x = variable, y = valueD, width = as.numeric(width), fill = donor, group = valueD, color = donor), 
             position = "fill", stat = "identity", show.legend = F) + 
    scale_fill_manual(name="", values=colors)+
    scale_color_manual(name="", values=colors)+
    geom_point(subset(Merge_DT, value != 0) , mapping = aes(x=variable, y=-0.2, color=donor, size = 60), size = 20)+ 
    #, guide = 'none'
    #scale_y_continuous(breaks=c(0, 0.5, 1), labels = percent, expand = c(0, 0)) + 
    scale_x_continuous(labels = c("0E0"="0", "2E6"="2", "4E6"="4", "6E6"="6", "8E6"="8")) + 
    #scale_x_continuous(labels = function(x) format(x, scientific = FALSE))+
    coord_cartesian(xlim = c(0,8206800), ylim = c(-0.3,1)) + 
    labs(title = " ", y = " ", x= "Chromosome position (Mb)\n",  colour = "")+
    facet_grid(Chr ~ ., scales = "free", space = "free")+
    #facet_matrix(rows = vars(Chr, Chr))+
    theme(plot.title = element_text(hjust = 1.5), 
          legend.position = "bottom", 
          legend.title=element_text(face="bold", size=30), 
          legend.text=element_text(face="bold", size=60),
          legend.direction = "horizontal",
          legend.key = element_rect(fill = "white"),
          axis.title.x = element_text(face="bold", size = 90, margin = margin(1, 21, 21, 21), vjust = -2),
          axis.text.x=element_text(face="bold", angle=0,  size = 60, colour = "black"), 
          axis.ticks.length.x = unit(1,"cm"),
          axis.ticks.x = element_line(size = 1, colour = "black"),
          #axis.text.x=element_text(size=rel(1.4), face="bold", colour = "black")
          axis.text.y = element_blank(),
          # axis.text.y= element_text(face="bold", angle=0,  size = 30, colour = "black"),
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(size=80, face = "bold", margin = margin(0,1,0,1, "cm")),
          strip.background=element_rect(fill = "grey", colour = "black", size = 0.7),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.7)        
    ) + 
    guides(colour=guide_legend(nrow=1)) 

)

dev.off()









