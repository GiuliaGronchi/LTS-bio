#!/usr/bin/Rscript

txt_to_r <- function(pos.txt) {
  pos <- read.table(pos.txt, dec = ",", header = FALSE)
  colnames(pos) <- c("lambda", "fondo_a", "valdec_a", "errdec_a", "fondo_b", "valdec_b", "errdec_b")
  pos$valore_a <- pos$valdec_a*pos$fondo_a/10
  pos$errore_a <- pos$errdec_a*pos$fondo_a/10
  pos$valore_b <- pos$valdec_b*pos$fondo_b/10
  pos$errore_b <- pos$errdec_b*pos$fondo_b/10
  pos$valore_r <- pos$valore_a/pos$valore_b #r è il rapporto tra canale A e B, a prescindere da posizione
  pos$errore_r <- pos$valore_r*sqrt((pos$errore_a/pos$valore_a)^2 + (pos$errore_b/pos$valore_b)^2)
  pos$errlog_r <- pos$valore_r*abs((pos$errore_a/pos$valore_a) - (pos$errore_b/pos$valore_b))
  return(pos)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Servono due argomenti: <file_pos1>.txt <file_pos2>.txt\n")
}

pos1 <- txt_to_r(args[1])
pos3 <- txt_to_r(args[2])

valore_tras <- sqrt(pos1$valore_r/pos3$valore_r)
errore_tras <- 0.5*valore_tras*sqrt((pos1$errore_r/pos1$valore_r)^2 + (pos3$errore_r/pos3$valore_r)^2)
errlog_tras <- 0.5*valore_tras*abs((pos1$errore_r/pos1$valore_r) - (pos3$errore_r/pos3$valore_r))

valore_ext <- -log(valore_tras)/10^7
errore_ext <- errore_tras/valore_tras/10^7
errlog_ext <- errlog_tras/valore_tras/10^7 #in questo caso assurdo usare formula angelo perché ext è funzione di una sola variabile

tab_ext <- cbind(pos1$lambda, valore_ext, errore_ext)
tablog_ext <- cbind(pos1$lambda, valore_ext, errlog_ext)
parti <- strsplit(args[1], "_pos", fixed = T)
write.table(tab_ext, paste("EXT_", parti[[1]][1], ".dat", sep = ""), sep = "\t", row.names = F, col.names = F)
write.table(tablog_ext, paste("EXT_", parti[[1]][1], "_errlog.dat", sep = ""), sep = "\t", row.names = F, col.names = F)
