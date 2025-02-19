#!/usr/bin/Rscript

#testato il 2018-04-24

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Servono 7 argomenti: <radice>_<ns>_<smax>.<suff> ns1 ns2 nsincr smax1 smax2 smaxincr\n")
}

radice_suff <- sub("_ns[0-9a-z_]+", "", args[1], perl = T) #radice di "<radice>_<ns>_<smax>.<suff>" di cui calcolare le concentrazione
parti <- strsplit(radice_suff, ".", fixed = T)
radice <- parti[[1]][1]
suff <- parti[[1]][2]

param <- strtoi(args[2:7]) #conversione a numeric
ns1 <- param[1]
ns2 <- param[2]
nsincr <- param[3]
smax1 <- param[4]
smax2 <- param[5]
smaxincr <- param[6]

nsnum <- (ns2 - ns1)/nsincr + 1
smaxnum <- (smax2 - smax1)/smaxincr + 1
ns <- sort(rep(seq(ns1, ns2, nsincr), smaxnum)) #fisso ns e scorro tutti gli smax: stesso ordine doppio for seguente
smax <- rep(seq(smax1, smax2, smaxincr), nsnum) 
concentrazione <- double(nsnum*smaxnum) #alloco prima per velocizzare parecchio il for

for (i in seq(ns1, ns2, nsincr)) {
	for (j in seq(smax1, smax2, smaxincr)) {
		dist <- read.table(paste(radice, "_ns", i, "_smax", j, ".", suff, sep = ""), header = FALSE)

		bin <- dist$V1[2] - 1 #le distribuzioni cominciano sempre da 1 e non da 0
		
		#i = ns1 + nsincr*(k-1), dove k è il vero indice che parte da 1 che voglio per concentrazione; analogo per j e smax
		k <- (i-ns1)/nsincr + 1
		l <- (j-smax1)/smaxincr + 1
		concentrazione[(k-1)*smaxnum+l] <- sum(bin*dist$V2)*10^21
	}
}

tab_conc <- cbind(ns, smax, concentrazione)
write.table(tab_conc, paste(sub("[A-Za-z]+", "CONC", radice, perl = T), ".dat", sep = ""), sep = "\t", row.names = F, col.names = F)

