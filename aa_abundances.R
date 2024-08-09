# ヴァネッサ
# (\___/)
# (='.'=)
# (")_(")

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


#---

aBunDances <- function(file){
      
      j <- 1
      aac <- data.frame(matrix(ncol=5, nrow=0))
      colnames(aac) <- c("Protein", "AA", "count", "pc", "ratio")
      for (i in file){
            data_1 <- aa_abundance(file[[j]])
            trans <- transmute(data_1, ratio = (data_1$pc / abundance_table$pc))
            data_1 <- cbind(data_1, trans)
            data_1$Protein <- rep(names(file[j]), times = nrow(data_1))
            data_1 <- data_1[,c(5, 1, 2, 3, 4)]
            aac <- rbind(aac, data_1)
            
            j <- j+1
      }
      
      res <- aac %>% group_by(AA) %>%
            summarise_at(vars(ratio), list(mean_rel_freq = mean))
      res %<>% as.data.frame()
      return(res)
      
}