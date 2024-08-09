#--simplified fasta reader

reader <- function (file) {

            lines <- scan(file, what = "", sep = "\n", quiet = TRUE)
            comments <- grep("^;", lines)
            
                  if (length(comments) > 0) {
                        lines <- lines[-comments]
            }
            ind <- which(substr(lines, 1L, 1L) == ">") # search lines with > as mark of new id
            nseq <- length(ind)
            if (nseq == 0) {
                  stop("incorrect file format")
            }
            start <- ind + 1
            end <- ind - 1
            end <- c(end[-1], length(lines))
            sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                                 collapse = ""))

            nomseq <- lapply(seq_len(nseq), function(i) {
                        firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
                        substr(firstword, 2, nchar(firstword))
                  })

            

            sequences <- lapply(sequences, s2c)
            
            
            names(sequences) <- nomseq
            return(sequences)
}