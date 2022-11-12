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

# grab a dataset from the ShinyHaplotypes output folder (too large to share through standard data repository)

df <- read.table(paste("~/SLIDES/", "chr1", ".", "ATCC64557", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

strains <- row.names(df)

strains[[length(strains)+1]] <- "ATCC64557"

sstrains <- sort(strains)

plotLength <- 3500000

vec <- list.files("~/SLIDES/")

vec <- gsub("\\..*", "", vec)

chromosomes <- unique(vec)

plotLength <- 600

# read a list of clade names

cladeName <- read.table("~/Downloads/cladeNames.csv", header = FALSE)

cladeNames <- as.character(cladeName[,1])

lineages <- c("B1", "B2", "B3", "B4", "C", "Ec",  "E1", "E2", "E3", "Er", "Ha", "L", 
              "L2", "L3", "Le", "Lee", "Lu", "M", "O", "P", "Pan", "S", "St", "T", "U")


# specify clade colors using old lineage identifiers (B was Brachiaria, now Urochloa; P was Paspalum, now Panicum)

#colors <- c(B1= "#0000FF", B2="#4169E1", B3= "#7B68EE", B4= "#87CEFA", C="#00FF00", Ec="#669999", 
#            E1= "#bf80ff", E2="#9999ff", E3="#CC99FF", Er="#CC6666", Ha = "#32cc99", 
#            L = "#a64dff", L2= "#619CFF", L3="#666600", Le="#9966CC", Lee="#9932CC", Lu= "#CC3333", M= "#cc9932", O = "#008000",
#           P="#FF6600", Pan="#336600", S= "#00BFC4", St= "#FF00FF", T= "#CC0000", U = "#F8766D")

myPlots <- list()


## code to replace html input

ui <- fluidPage(
  #  tags$head(tags$style(HTML(".titlePanel {font-size: 5px;}"))),
  #  titlePanel("Select strain, clade(s) and chromosome", windowTitle = "Haplotype Divergence"),
  title = "Haplotype Divergence",
  
  plotOutput("Windows"),
  
  hr(),
  
  fluidRow(
    
    column(5,
           radioButtons(inputId = "inputStrain", "Select strain:", 
                        inline = TRUE, choices = sstrains,
                        selected = sstrains[[1]]),
           br(),
           sliderInput("yaxis", "Haplotype divergence:",
                       min = 0, max = 1,
                       step=0.1, value = c(0,0.5)),
           br(),
           actionButton("print", "PrintPlot")
    ),
    
    column(6, offset = 1,
           selectInput("comparison", "Type of comparison",
                       choices = c("against lineages", "against individual strains")),
           conditionalPanel(
             condition = "input.comparison == 'against lineages'",
             checkboxGroupInput("inputLineage", "Select lineages(s) for comparison:", 
                                inline = TRUE, choices = lineages, # cladeNames,
                                selected = 'B1') #cladeNames[[1]]),
           ),
           conditionalPanel(
             condition = "input.comparison == 'against individual strains'",
             checkboxGroupInput("compareStrain", "Select strains(s) for comparison:", 
                              inline = TRUE, choices = strains, # cladeNames,
                              selected = strains[[2]]) #cladeNames[[1]]),
           )
    )
  )
)


server <- function(input, output, session) {
  
  region <- reactive(  {input$region}   )
  output$Windows <- renderPlot({
    # Read in haplotype information
    
  observeEvent(input$print, {
    session$sendCustomMessage(type = 'testmessage',
                              message = 'Thank you for clicking')
  })
      
    colors <- c(T="#000000", B1="#D55E00", O="#0072B2", S= "#0072B2", E1="#009E73", S="#F0E442", Pan="#0072B2", H= "#999999", B3= "#AA4499")
    
    for (f in 1:7) {
      
#    df <- read.table(paste("~/SLIDES/",  "chr", f, ".", input$inputStrain, ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)  
    df <- read.table(paste("~/SLIDES/",  "chr", f, ".", "PY0925", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)  
    
    chrSize= length(colnames(df))-1
    # Get sliding window information

#    df2 <- df %>% rownames_to_column('Strains') %>%  {if(input$comparison == 'against individual strains') filter(., Strains %in% input$compareStrain) else filter(., clade %in% input$inputLineage)} %>% column_to_rownames('Strains') # c("B3", "E2")
    df2 <- df %>% rownames_to_column('Strains') %>% filter(., Strains %in% c("Arcadia2", "Guy11", "B51", "Bm88324", "U169-v1", "Br35", "Pr8202", "FPH-2015-44")) %>% column_to_rownames('Strains') # c("B3", "E2")
    end <- ncol(df2)
    
    win_pos <- c(1:end)
    colnames(df2) <- win_pos
    df2 <- df2 %>% rownames_to_column('Strains') 
    colnames(df2)[end+1] <- "clade"
    
    df3<- df2 %>% select(Strains, clade, win_pos)%>%gather(key = "win_pos", value = "divergence", -Strains, -clade) 
    
    if (f == 7) {
      
      myPlots[[f]] <- ggplot(data = df3, aes(x = as.numeric(win_pos), y = as.numeric(divergence))) + 
        ggtitle(paste("Chr", f, sep="")) +
        theme(plot.title=element_text(hjust=0.98, vjust=0.98, face='bold', size = 12, margin=margin(t=20,b=-20))) + 
        geom_line(aes(color = clade, group = Strains), size = 0.8) +
        xlim(0,1200) + ylim(0,0.4) +
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
      geom_line(aes(color = clade, group = Strains), size = 0.8)+
      #labs(title = "Haplotype Divergence", y = "Divergence (%)", x= "Window position")+
      xlim(0,1100) + 
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
    
    pdf("B71_ShinyHaplotypes.pdf", 8, 16)
    p<-grid.arrange(grobs = myPlots, ncol = 1, nrow = 7, heights = c(1,1,1,1,1,1,1.2), 
      #top = textGrob("Haplotype Divergence: B71 x B51+Bm88324+Guy11+Pr8202+U169+Up35", gp = gpar(fontface = "bold", cex = 1.2)),
      left = textGrob("Divergence", rot = 90),
      bottom = textGrob("Window Position"))
    print(p)
    dev.off()
    
  })
  
}

shinyApp(ui = ui, server = server)






