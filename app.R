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

# define user interface

ui <- page_sidebar(
      title = "aBunDances",
      sidebar = sidebar("check this out!", position = "right"),
      "Calculate relative amino acid composition",
      card(
            card_header("Getting started"),
            "To get started, download the fasta files of your POIs at
            https://www.uniprot.org/id-mapping",
            card_image("bun.jpeg", width = "500px", height = "300px")
      ),
   
      value_box(
            title = "boxV",
            value = 100,
            showcase = bsicons::bs_icon("bar-chart"),
            theme = "teal"
      ),
      card(
            card_header("Upload fasta file"),
            fileInput(inputId = "file", label = NULL, accept = ".fasta")
      ),
      
      card(
            actionButton("calculate", label = "Calculate!",
                         icon = icon("jedi"))
      ),
      
      card(
            
            dataTableOutput("results")
            
      )
      
)



# define server logic
server <- function(input, output) {
     
      
      aa_abundance <- function(x){
            
            A_Ala <- sum(str_count(x, "M"))
            C_Cys <- sum(str_count(x, "C"))
            D_Asp <- sum(str_count(x, "D"))
            E_Glu <- sum(str_count(x, "E"))
            F_Phe <- sum(str_count(x, "F"))
            G_Gly <- sum(str_count(x, "G"))
            H_His <- sum(str_count(x, "H"))
            I_Ile <- sum(str_count(x, "I"))
            K_Lys <- sum(str_count(x, "K"))
            L_Leu <- sum(str_count(x, "L"))
            M_Met <- sum(str_count(x, "M"))
            N_Asn <- sum(str_count(x, "N"))
            P_Pro <- sum(str_count(x, "P"))
            Q_Gln <- sum(str_count(x, "Q"))
            R_Arg <- sum(str_count(x, "R"))
            S_Ser <- sum(str_count(x, "S"))
            T_Thr <- sum(str_count(x, "T"))
            U_Sec <- sum(str_count(x, "U"))
            V_Val <- sum(str_count(x, "V"))
            W_Trp <- sum(str_count(x, "W"))
            Y_Tyr <- sum(str_count(x, "Y"))
            
            
            ab_table <- data.table::transpose(data.frame(
                  A_Ala, C_Cys, D_Asp, E_Glu, F_Phe, G_Gly, H_His, I_Ile, K_Lys,
                  L_Leu, M_Met, N_Asn, P_Pro, Q_Gln, R_Arg, S_Ser, T_Thr, U_Sec,
                  V_Val, W_Trp, Y_Tyr), keep.names = "AA")
            ab_table <- rename(ab_table, "count" = "V1")
            total_count <- sum(ab_table$count)
            freq <- transmute(ab_table, pc = (count / total_count)*100)
            ab_table <- cbind(ab_table, freq)
            
            return(ab_table)
      }
      
      output$results <- renderTable({ observeEvent(input$calculate, {
            req(input$file)
            #file <- input$file
                  ext <- tools::file_ext(file$datapath)
                  #req(file)
                  validate(need(ext == "fasta", "Please upload fasta file"))
                  
                  data <- read.fasta(file$datapath, seqtype = "AA", whole.header = FALSE,
                                     set.attributes = FALSE)
                  abundance_table <- read.table("AA_frequencies.csv")
                  
                  
                  j <- 1
                  aac <- data.frame(matrix(ncol=5, nrow=0))
                  colnames(aac) <- c("Protein", "AA", "count", "pc", "ratio")
                  for (i in data){
                        data_1 <- aa_abundance(data[[j]])
                        trans <- transmute(data_1, ratio = (data_1$pc / abundance_table$pc))
                        data_1 <- cbind(data_1, trans)
                        data_1$Protein <- rep(names(data[j]), times = nrow(data_1))
                        data_1 <- data_1[,c(5, 1, 2, 3, 4)]
                        aac <- rbind(aac, data_1)
                        
                        j <- j+1
                  }
                  
                  res <- aac %>% group_by(AA) %>%
                        summarise_at(vars(ratio), list(mean_rel_freq = mean))
                  res %<>% as.data.frame()
                  return(res)
                  
            })
      })
      

}


# run the app
shinyApp(ui = ui, server = server)








# &c-7BvMjGeWcM@Q


# output$plot <- renderPlot({
#       
#       pl <- ggplot(output$results, aes(x = AA, y = mean_rel_freq))
#       pl1 <- pl + geom_col(aes(), position = "dodge", color = "black",
#                            width = 0.7, fill = "darkred") + ylim(0,1.6)
#       pl2 <- pl1 + theme_clean() + ggtitle("Amino acid composition",
#                                            subtitle = "decreased measurements - mean values") + 
#             ylab("Relative Frequency") + theme(
#                   axis.title.x = element_text(size = 12, face = "bold"),
#                   axis.title.y = element_text(size = 12, face = "bold"),
#                   axis.text.x = element_text(face = "bold"),
#                   axis.text.y = element_text(size = 10, face = "bold")
#             ) + 
#             scale_fill_manual(values = c("red"))
#       pl3 <- pl2 + geom_hline(yintercept = 1, linetype = "dashed", size = 1.2)
#       plot <- pl3 + geom_text(aes(label = round(mean_rel_freq,2)), vjust = -0.2)
#       plot
#       
# })
# 



#-------------------------------------------------------------------------------

## human proteom fasta file read out
# calculate relative amino acid abundance

library(seqinr)
library(stringr)
library(data.table)
library(tidyverse)
library(ggthemes)
library(RColorBrewer)
library(data.table)

# creating reference file
###########
# fasta downloaded 12.10.2022, UniProt DB, Swiss-Prot curated entries only
file <- read.fasta(
      "C:/Users/Beutgen/Documents/R_projects/beta_project/data/human_proteome_fasta/human_proteome.fasta",
      seqtype = "AA" ,whole.header = FALSE, set.attributes = FALSE)

A_Ala <- sum(str_count(file, "M"))
C_Cys <- sum(str_count(file, "C"))
D_Asp <- sum(str_count(file, "D"))
E_Glu <- sum(str_count(file, "E"))
F_Phe <- sum(str_count(file, "F"))
G_Gly <- sum(str_count(file, "G"))
H_His <- sum(str_count(file, "H"))
I_Ile <- sum(str_count(file, "I"))
K_Lys <- sum(str_count(file, "K"))
L_Leu <- sum(str_count(file, "L"))
M_Met <- sum(str_count(file, "M"))
N_Asn <- sum(str_count(file, "N"))
P_Pro <- sum(str_count(file, "P"))
Q_Gln <- sum(str_count(file, "Q"))
R_Arg <- sum(str_count(file, "R"))
S_Ser <- sum(str_count(file, "S"))
T_Thr <- sum(str_count(file, "T"))
U_Sec <- sum(str_count(file, "U"))
V_Val <- sum(str_count(file, "V"))
W_Trp <- sum(str_count(file, "W"))
Y_Tyr <- sum(str_count(file, "Y"))

total_count <- sum(A_Ala, C_Cys, D_Asp, E_Glu, F_Phe, G_Gly, H_His, I_Ile,
                   K_Lys, L_Leu, M_Met, N_Asn, P_Pro, Q_Gln, R_Arg, S_Ser,
                   T_Thr, U_Sec, V_Val, W_Trp, Y_Tyr)

abundance_table <- data.table::transpose(data.frame(
      A_Ala, C_Cys, D_Asp, E_Glu, F_Phe, G_Gly, H_His, I_Ile, K_Lys,
      L_Leu, M_Met, N_Asn, P_Pro, Q_Gln, R_Arg, S_Ser, T_Thr, U_Sec,
      V_Val, W_Trp, Y_Tyr), keep.names = "AA")
abundance_table <- rename(abundance_table, "count" = "V1")
rel_freq <- transmute(abundance_table, pc = (count / total_count)*100)
abundance_table <- cbind(abundance_table, rel_freq)
abundance_table$pc <- round(abundance_table$pc, 3)

write.table(abundance_table,
            file = "C:/Users/Beutgen/Documents/R_projects/a_bun_dances/AA_frequencies.csv")
###########
# only soma targets as background
abundance_table <- read.table("C:/Users/Beutgen/Documents/R_projects/beta_project/data/AA_freq_soma_targets.csv")

# whole human proteome as background
abundance_table <- read.table("C:/Users/Beutgen/Documents/R_projects/beta_project/data/AA_frequencies.csv")
abundance_table$pc <- round(abundance_table$pc, 3)

aa_abundance <- function(x){
      
      A_Ala <- sum(str_count(x, "M"))
      C_Cys <- sum(str_count(x, "C"))
      D_Asp <- sum(str_count(x, "D"))
      E_Glu <- sum(str_count(x, "E"))
      F_Phe <- sum(str_count(x, "F"))
      G_Gly <- sum(str_count(x, "G"))
      H_His <- sum(str_count(x, "H"))
      I_Ile <- sum(str_count(x, "I"))
      K_Lys <- sum(str_count(x, "K"))
      L_Leu <- sum(str_count(x, "L"))
      M_Met <- sum(str_count(x, "M"))
      N_Asn <- sum(str_count(x, "N"))
      P_Pro <- sum(str_count(x, "P"))
      Q_Gln <- sum(str_count(x, "Q"))
      R_Arg <- sum(str_count(x, "R"))
      S_Ser <- sum(str_count(x, "S"))
      T_Thr <- sum(str_count(x, "T"))
      U_Sec <- sum(str_count(x, "U"))
      V_Val <- sum(str_count(x, "V"))
      W_Trp <- sum(str_count(x, "W"))
      Y_Tyr <- sum(str_count(x, "Y"))
      
      
      ab_table <- data.table::transpose(data.frame(
            A_Ala, C_Cys, D_Asp, E_Glu, F_Phe, G_Gly, H_His, I_Ile, K_Lys,
            L_Leu, M_Met, N_Asn, P_Pro, Q_Gln, R_Arg, S_Ser, T_Thr, U_Sec,
            V_Val, W_Trp, Y_Tyr), keep.names = "AA")
      ab_table <- rename(ab_table, "count" = "V1")
      total_count <- sum(ab_table$count)
      freq <- transmute(ab_table, pc = (count / total_count)*100)
      ab_table <- cbind(ab_table, freq)
      
      return(ab_table)
}


## calculate fold changes in relative frequency
# down = significantly decreased measurements with FDR < 0.05
down <- read.fasta("C:/Users/Beutgen/Documents/R_projects/beta_project/data/newLists/down_list.fasta",
                   seqtype = "AA", whole.header = FALSE, set.attributes = FALSE)


j <- 1
aac_decreased <- data.frame(matrix(ncol=5, nrow=0))
colnames(aac_decreased) <- c("Protein", "AA", "count", "pc", "ratio")
for (i in down){
      down_1 <- aa_abundance(down[[j]])
      xyz <- transmute(down_1, ratio = (down_1$pc / abundance_table$pc))
      down_1 <- cbind(down_1, xyz)
      down_1$Protein <- rep(names(down[j]), times = nrow(down_1))
      down_1 <- down_1[,c(5, 1, 2, 3, 4)]
      aac_decreased <- rbind(aac_decreased, down_1)
      
      j <- j+1
}

down_mean <- aac_decreased %>%
      group_by(AA) %>%
      summarise_at(vars(ratio), list(mean_rel_freq = mean))


# plots
# decreased single
pl <- ggplot(aac_decreased, aes(x = AA, y = ratio))
pl1 <- pl + geom_col(aes(fill=Protein), position = "dodge", color = "black",
                     width = 0.7) + ylim(0,3)
pl2 <- pl1 + theme_clean() + ggtitle("Amino acid composition",
                                     subtitle = "decreased measurements") + 
      ylab("Relative Frequency") + 
      scale_fill_brewer(type = "seq", palette = "Spectral", direction = 2,
                        aesthetics = "fill") + theme(
                              axis.title.x = element_text(size = 12, face = "bold"),
                              axis.title.y = element_text(size = 12, face = "bold"),
                              axis.text.x = element_text(face = "bold"),
                              axis.text.y = element_text(size = 10, face = "bold")
                        )
pl2 + geom_hline(yintercept = 1, linetype = "dashed", size = 1.2)

# decreased mean
pl <- ggplot(down_mean, aes(x = AA, y = mean_rel_freq))
pl1 <- pl + geom_col(aes(), position = "dodge", color = "black",
                     width = 0.7, fill = "darkred") + ylim(0,1.6)
pl2 <- pl1 + theme_clean() + ggtitle("Amino acid composition",
                                     subtitle = "decreased measurements - mean values") + 
      ylab("Relative Frequency") + theme(
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(face = "bold"),
            axis.text.y = element_text(size = 10, face = "bold")
      ) + 
      scale_fill_manual(values = c("red"))
pl3 <- pl2 + geom_hline(yintercept = 1, linetype = "dashed", size = 1.2)
pl3 + geom_text(aes(label = round(mean_rel_freq,2)), vjust = -0.2)

