#!/usr/bin/Rscript

#calcola la differenza percentuale tra due colonne di dati (tipo errori calcolati in modo diverso)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Servono due argomenti: <file_dat1> <file_dat2>\n")
}

dat1 <- read.table(args[1])
dat2 <- read.table(args[2])
colnames(dat1) <- c("x", "y", "err_y")
colnames(dat2) <- c("x", "y", "err_y")
diff <- abs(dat1$err_y - dat2$err_y)/dat1$err_y*100 #percentuale
table <- cbind(dat1$x, diff)
write.table(table, "diff.txt", sep = "\t", row.names = F, col.names = F)
