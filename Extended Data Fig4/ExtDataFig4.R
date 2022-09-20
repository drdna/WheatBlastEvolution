# grab dependencies

'%!in%' <- function(x,y)!('%in%'(x,y))
require('ggplot2')
library(data.table)
library(shiny)
library(dplyr)
library(RColorBrewer)
library("tidyverse")
library(tibble)  # for `rownames_to_column` and `column_to_rownames`
library(magrittr)
library(grid)
library(gridExtra)
#library("pandas")

# read in color set

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# grab initial dataset

df <- read.table(paste("~/NEE_SHINY/", "Chr1", ".", "ATCC64557", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
df <- df[-ncol(df)]


# grab a list of strain IDs

strains <- row.names(df)

strains[[length(strains)+1]] <- "ATCC64557"

sstrains <- sort(strains)


# grab a list of chromosome IDs

vec <- list.files("~/NEE_SHINY/")

vec <- gsub("\\..*", "", vec)

chromosomes <- unique(vec)


# read a list of clade names

cladeName <- read.table("~/strain.idfile", header = FALSE)

cladeNames <- as.character(cladeName[,1])

Lineages <- c("B1", "B2", "B3", "B4", "C", "Ec",  "E1", "E2", "E3", "Er", "Ha", "L", 
              "L2", "L3", "Le", "Lee", "Lu", "M", "O", "P", "Pan", "S", "St", "T", "U")


# specify clade colors

colors <- c(U1= "#D55E00", U2="#0072B2", U3= "#AA4499", C="black", Ec="black", 
            E1="#009E73", E2="#009E73", E3="#009E73", Er="black", H = "#FF00FF", L = "#a64dff", L2= "black", 
            Le="#0072B2", Lee="#0072B2", Lu= "#999999", M= "black", 
            O="#0072B2", P="black", P2="#0072B2", 
            S= "#0072B2", St= "#AA4499", T="#0000FF", U = "#AA4499")


# empty out list of plots (for creating new plots each time script is run)

myPlots <- list()


# generate plots for chromosomes 1 through 7

for (f in 1:7) {
  
  
  # modify next line to select a different reference strain  
  df <- read.table(paste("~/NEE_SHINY/",  "Chr", f, ".", "B71", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)  
  df <- df[-ncol(df)]     # strip off last column - redundant and messes up color assignment
  
  chrSize= length(colnames(df))-1
  # Get sliding window information
  
  # grab rownames and put in a new column named strains - filter by strains
  
  
  # uncomment/modify next line to select specific strain comparisons
  #df2 <- df %>% rownames_to_column('Strains') %>% filter(., Strains %in% c("ATCC64577", "CHRF", "Ds555i", "PtKY18-1", "PL3-1", "TF05-1", "PY86", "pg1213-22", "P25", "P28", "U234", "T47-3", "T42-2", "T12-8", "T2-1", "T21-1", "T13-3", "T37-2", "T50-3", "Br7", "Br80", "Py221", "T46-2", "PY6017", "BR118", "Br130", "P3", "PY5033", "PY36", "PY0925", "WB032i" ))
  df2 <- df %>% rownames_to_column('Strains') %>% filter(., clade %in% c("U1", "U2", "U3", "E1", "E2", "E3", "Er", "H", "L2", "Le", "Lee", "Lu", "M", "O", "P", "P2", "S", "St", "U"))
  # uncomment/modify next line for specific lineage comparisons
  #df2 <- df %>% rownames_to_column('Strains') %>% filter(., clade %in% c("L", "T"))
  end <- ncol(df2)
  
  
  # define position of window
  win_pos <- c(1:end)
  colnames(df2) <- win_pos
  df2 <- df2 %>% rownames_to_column('Strains') 
  colnames(df2)[end+1] <- "clade"
  
  
  # gather data for plotting
  df3<- df2 %>% select(Strains, clade, win_pos)%>%gather(key = "win_pos", value = "divergence", -Strains, -clade) 
  
  if (f == 7) {
    
    myPlots[[f]] <- ggplot(data = df3, aes(x = as.numeric(win_pos), y = as.numeric(divergence))) + 
      ggtitle(paste("Chr", f, sep="")) +
      theme(plot.title=element_text(hjust=0.98, vjust=0.98, face='bold', size = 12, margin=margin(t=20,b=-20))) + 
      geom_line(aes(color = clade, group = Strains), size = 0.2, linetype = "dotted") +
      ylim(0,0.4) +
      scale_x_continuous(limits=c(0, 2250), breaks = seq(from = 0, to = 2250, by = 100)) +
      scale_color_manual(values = colors) +
      theme(plot.margin = margin(0.5,1,1.5,5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 0.3),
            axis.text.x=element_text(size=rel(0.8), face="bold", colour = "black"),
            axis.text.y=element_text(size=rel(0.8), face="bold", colour = "black"),
            axis.title=element_blank(),
            axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"),
            legend.title=element_blank(),
            legend.position = "none")
  }
  
  else {
    
    myPlots[[f]] <- ggplot(data = df3, aes(x = as.numeric(win_pos), y = as.numeric(divergence))) + 
      ggtitle(paste("Chr", f, sep="")) +
      theme(plot.title=element_text(hjust=0.98, vjust=0.98, face='bold', size = 12, margin=margin(t=20,b=-20))) + 
      geom_line(aes(color = clade, group = Strains), size = 0.2, linetype = "dotted")+
      #labs(title = "Haplotype Divergence", y = "Divergence (SNPs/200 variant sites)", x= "Window position")+
      scale_x_continuous(limits=c(0, 2250), breaks = seq(from = 0, to = 2250, by = 100)) +
      ylim(0,0.4) +
      scale_color_manual(values= colors) +
      theme(plot.margin=margin(0.5,1,1.5,5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 0.3),
            axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=rel(0.8), face="bold", colour = "black"),
            axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"),
            legend.title=element_blank(),
            legend.position = "none")
  }
  
}


# uncomment next section to print to file
pdf("B71vAll.pdf", 8.5, 11)
p<-grid.arrange(grobs = myPlots, ncol = 1, nrow = 7, heights = c(1,1,1,1,1,1,1.2), 
                top = textGrob("Haplotype Divergence: B71 versus All PoT/PoL1", gp = gpar(fontface = "bold", cex = 1.2)),
                left = textGrob("Divergence (SNPs/200 variant sites)", rot = 90),
                bottom = textGrob("Window Position"))
print(p)
dev.off()









