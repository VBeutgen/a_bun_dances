# ヴァネッサ
# (\___/)
# (='.'=)
# (")_(")

# Setup

# Source scripts, load libraries, and read data sets at the beginning of app.R
# outside of the server function. Shiny will only run this code once, which is
# all you need to set your server up to run the R expressions contained in server.

library(shiny)
library(bslib)
library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)
library(ggthemes)
library(RColorBrewer)
library(data.table)
library(magrittr)

source("aa_abundances.R")
source("reader.R")

abundance_table <- read.table("data/AA_frequencies.csv")


# define user interface

ui <- page_sidebar(
      theme = bs_theme(preset = "vapor"),
      title = "aBunDances",
      sidebar = sidebar(
            "Data Import / Export",
            position = "right",
            fileInput("file", label = "Upload fasta file", accept = ".fasta"),
            actionButton("calculate", label = "Calculate!", icon = icon("jedi"))
            ),
      "Calculate relative amino acid composition",
      card(
            card_header("Getting started"),
            "To get started, download the fasta files of your POIs at
            https://www.uniprot.org/id-mapping",
            card_image("data/bun.jpeg", width = "500px", height = "300px")
      ),
   
      value_box(
            title = "showing results for",
            value = textOutput("name"),
            showcase = bsicons::bs_icon("bar-chart")
      ),
      
      card(
            plotOutput("plot")
      ),
      
      card(
            tableOutput("results")
      )
      
)

server <- function(input, output) {
      
      # Reactive file input
      dataInput <- reactive({
            req(input$file)
            
            
            # Normalize the file path
            filepath <- normalizePath(input$file$datapath, winslash = "/")
            print(filepath)
            
            # Check if the file exists
            if (!is.null(filepath) && file.exists(filepath)) {
                  print("File exists and is accessible.")
                  
                  # Read the fasta file
                  fasta_data <- reader(filepath)
                  return(fasta_data)
            } else {
                  print("File does not exist or path is NULL")
                  return(NULL)
            }
      })
      
      abun <- eventReactive(input$calculate, {
            req(dataInput())
            aBunDances(dataInput())
      })

      output$results <- renderTable({
            abun()
      })
      
      output$plot <- renderPlot({
            req(abun())
            
            pl <- ggplot(abun(), aes(x = AA, y = mean_rel_freq))
            pl1 <- pl + geom_col(aes(), position = "dodge", color = "black",
                                 width = 0.7, fill = "lightgreen") + ylim(0,1.6)
            pl2 <- pl1 + theme_clean() + ggtitle("Amino acid composition") + 
                  ylab("Relative Frequency") + theme(
                        axis.title.x = element_text(size = 12, face = "bold"),
                        axis.title.y = element_text(size = 12, face = "bold"),
                        axis.text.x = element_text(face = "bold"),
                        axis.text.y = element_text(size = 10, face = "bold")
                  ) + 
                  scale_fill_manual(values = c("red"))
            pl3 <- pl2 + geom_hline(yintercept = 1, linetype = "dashed", size = 1.2)
            pl3 + geom_text(aes(label = round(mean_rel_freq,2)), vjust = -0.2)
      })
      
      output$name <- 
            renderPrint({
                  input$file$name
            })
}

# Run the app
shinyApp(ui = ui, server = server)