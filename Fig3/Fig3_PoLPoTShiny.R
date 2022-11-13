# grab dependencies

'%!in%' <- function(x,y)!('%in%'(x,y))
require('ggplot2')
library(data.table)
library(knitr)
library(shiny)
library(dplyr)
library(RColorBrewer)
library("tidyverse")
library(tibble)  # for `rownames_to_column` and `column_to_rownames`
library(magrittr)
#library("pandas")


# grab clade affiliation details

cladeNames <- read.table("~/NEE_CP/strain.idfile", header = FALSE)

cladeIDs <- as.character(cladeNames[,2])

lineages <- sort(unique(cladeIDs))

strainIDs <- cladeNames[,1]

strainAffil <- paste(as.character(cladeNames[,1]), " ", "(", cladeIDs, ")", sep = "")

names(strainIDs) <- strainAffil


# grab initial dataset

df <- read.table(paste("~/NEE_SHINY/", "Chr1", ".", "ATCC64557", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)

df <- df[-ncol(df)]

strains <- row.names(df)

strains[[length(strains)+1]] <- "ATCC64557"

sstrains <- sort(strains)

plotLength <- 2000

vec <- list.files("~/NEE_SHINY/")

vec <- gsub("\\..*", "", vec)

chromosomes <- unique(vec)

plotLength <- 600

# specify clade colors

colors <- c(U1= "#D55E00", U2="#4169E1", U3= "#FF00FF", Ho= "#87CEFA", C="#00FF00", Ec="#669999", 
            E1= "#009E73", E2="#9999ff", E3="#CC99FF", Er="#CC6666", H = "#32cc99", 
            L = "#a64dff", L2= "#619CFF", L3="#666600", Le="#9966CC", Lee="#9932CC", Lu= "#999999", M= "#cc9932", O= "#008000",
            P="#FF6600", P2="#336600", S= "#00BFC4", St= "#AA4499", "T"= "#0000FF", U = "#F8766D")


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
                        inline = TRUE, choices = strainIDs,
                        selected = strainIDs[[1]]),
           br(),
           sliderInput("yaxis", "Haplotype divergence:",
                       min = 0, max = 1,
                       step=0.1, value = c(0,0.5))
    ),
    
    column(6, offset = 1,
           selectInput("comparison", "Type of comparison",
                       choices = c("against lineages", "against individual strains")),
           conditionalPanel(
             condition = "input.comparison == 'against lineages'",
             checkboxGroupInput("inputLineage", "Select lineages(s) for comparison:", 
                                inline = TRUE, choices = lineages, # cladeNames,
                                selected = lineages[[1]]),
                                actionLink("selectall","Select All"),
           ),
           conditionalPanel(
             condition = "input.comparison == 'against individual strains'",
             checkboxGroupInput("compareStrain", "Select strains(s) for comparison:", 
                              inline = TRUE, choices = sstrains, # cladeNames,
                              selected = strains[[2]]) #cladeNames[[1]]),
           ),
           br(),
           radioButtons(inputId = "inputChromo", "Chromosome to plot:",
                        inline = TRUE, choices = chromosomes,
                        selected = chromosomes[[1]]),
           br(),
           sliderInput("region", "Select chromosome interval:",
                       min = 0, max = 2000,
                       step=1, value = c(1,2000))
    )
  )
)


server <- function(input, output, session) {
  
  region <- reactive(  {input$region}   )
  output$Windows <- renderPlot({
    # Read in haplotype information 
    df <- read.table(paste("~/NEE_SHINY/",input$inputChromo, ".", 
                           input$inputStrain, ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE) 
    df <- df[-ncol(df)]
    
    chrSize= length(colnames(df))-1
    # Get sliding window information
    df2 <- df %>% rownames_to_column('Strains') %>%  {if(input$comparison == 'against individual strains') filter(., Strains %in% input$compareStrain) else filter(., clade %in% input$inputLineage)} %>% column_to_rownames('Strains') # c("B3", "E2")
    end <- ncol(df2)
    win_pos <- c(1:end)
    colnames(df2) <- win_pos
    df2 <- df2 %>% rownames_to_column('Strains') 
    colnames(df2)[end+1] <- "clade"
  
    #updateSliderInput(session, "region", min = min(region), max = max(region))
    df3<- df2 %>% select(Strains, clade, win_pos)%>%gather(key = "variable", value = "value", -Strains, -clade) 
    updateSliderInput(session, inputId = "region",min=0, max = chrSize)
    ggplot(data = df3, aes(x = as.numeric(variable),y = as.numeric(value))) + 
      geom_line(aes(color = clade, group = Strains), size = 0.2)+
      labs(title = "Haplotype Divergence", y = "Divergence (%)", x= "Window position")+
      xlim(input$region)+ 
      ylim (input$yaxis)+
      scale_color_manual(values = colors) +
      guides(color = guide_legend(override.aes = list(size = 2))) + 
      theme(plot.title = element_text(size = 20, hjust = 0.5, face="bold"), panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black", size = 0.3),
            axis.text=element_text(size=rel(1.2), face="bold", colour = "black"),
            axis.line = element_line(colour = "black", size = 0.3, linetype = "solid"),
            axis.title = element_text(face="bold", size = 18, vjust = 3),
            legend.text= element_text(size=14, face="bold"),
            legend.title = element_text(colour="black", size=20, face="bold"),
            #legend.title=element_blank(),
            legend.direction = "horizontal",
            legend.position = "top",
            legend.justification = c(0.5, 1),
            legend.background = element_rect(fill="white",size=0.3, linetype="solid",colour ="black"),
            legend.key = element_rect(fill = "white")
            ) 
  })
}

shinyApp(ui = ui, server = server)






