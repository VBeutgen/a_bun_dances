# ヴァネッサ
# (\___/)
# (='.'=)
# (")_(")

# aBunDances

# Setup

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
library(shinycssloaders)

source("aa_abundances.R")
source("reader.R")




# define user interface

ui <- page_sidebar(
            theme = bs_theme(preset = "vapor"),
            title = "aBunDances",
            width = 3000,
            sidebar = sidebar("Data Import / Export",
                  
                  fileInput("file", label = "Upload fasta file",
                            accept = ".fasta"),
                  
                  actionButton("calculate", label = "Calculate!",
                               icon = icon("martini-glass-citrus"))
                  ),
      
      fluidPage(
            tags$h1("Calculate Relative Amino Acid Composition"),
            
            card(
                  card_header(tags$h2("Getting started")),
                  "To get started, download the FASTA file of your POIs at:",
                  tags$a("https://www.uniprot.org/id-mapping",
                  href = "https://www.uniprot.org/id-mapping")
            ),

            navset_card_underline(
                  
                  nav_panel("Plot",
                            withSpinner(
                              plotOutput("plot"),
                              type = 6, color = "purple"
                            )
            ),
            
                  nav_panel("Table",
                             tableOutput("results")
            ),
            
                  nav_panel("About",
                            
                            card(
                            " Human proteome reference data was downloaded in
                            FASTA format from UniProt (12.10.2022).
                            Swissprot-reviewed entries and canonical sequences
                            (no isoforms) were considered exclusively,
                            resulting in 20,398 protein sequences.
                            The retrieved sequence data was used to calculate
                            reference values for expected amino acid frequencies
                            normalized to 1 and are compared to the observed
                            frequencies of the uploaded protein list in the app.",
                            card_footer("Author:",
                                        tags$a("V.M. Beutgen",
                                               href = "https://github.com/VBeutgen")))
                            
                            )),
      
            
            card(
                  downloadButton("loadResults", "Download")
            )
      )
)


server <- function(input, output) {
      
      abundance_table <- read.table("data/AA_frequencies.csv")
      
      dataInput <- reactive({
            req(input$file)
            
            # Normalize the file path
            filepath <- normalizePath(input$file$datapath, winslash = "/")
            print(filepath)
            
            # Check if the file exists
            if (!is.null(filepath) && file.exists(filepath)) {
                  print("File exists and is accessible.")
                  
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
                                 width = 0.7, fill = "turquoise") + ylim(0,1.6)
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
      
      results <- eventReactive(input$calculate, {
            require(abun())
            abun()
      })

      output$loadResults <- downloadHandler(
            filename = function() {
            paste("results-", Sys.Date(), ".zip", sep = "")
      },
            content = function(zipfile) {
                  
                  temp_dir <- tempdir()
                  
                  csvfile <- file.path(temp_dir, "results.csv")
                  plotfile <- file.path(temp_dir, "plot.png")
                  
            write.csv(abun(), csvfile, row.names = F)
            
            ggsave(plotfile, plot = {
                  pl <- ggplot(abun(), aes(x = AA, y = mean_rel_freq))
                  pl1 <- pl + geom_col(aes(), position = "dodge", color = "black",
                                       width = 0.7, fill = "turquoise") + ylim(0,1.6)
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
            }, width = 10, height = 8)
            
            oldwd <- setwd(temp_dir)
            on.exit(setwd(oldwd))
            
            zip(zipfile, files = c("results.csv", "plot.png"))
      },
      contentType = "application/zip"
      )
      
}
      


# Run the app
shinyApp(ui = ui, server = server)