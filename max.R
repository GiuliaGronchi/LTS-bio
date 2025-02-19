#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Serve un argomento: LTS<param>.sol\n")
}

table <- read.table(args[1], header = FALSE)
max <- max(table$V2)
parti <- strsplit(args[1], "smax", fixed = T)
smax <- strsplit(parti[[1]][2], ".sol", fixed = T)
point <- c(strtoi(smax), max)
point
