#!/bin/bash

#2018-04-13, prima delle 9:30 - verificato che dà stesso output di par_poly_old.sh su poly_460nm_c10^9

#Permette di variare automaticamente ns (numero di punti) e smax (estremo superiore dell'intervallo dei raggi) in LTS.par

ns1=$1		#numero di punti iniziale
ns2=$2		#numero di punti finale
nsincr=$3	#incremento del numero di punti

smax1=$4	#estremo superiore iniziale
smax2=$5	#estremo superiore finale
smaxincr=$6	#incremento dell'estremo superiore

for ((i=$ns1; i<=$ns2; i+=$nsincr))
do
	sed -E -i "s/[0-9]+[ \t]+ns/$i\tns/" LTS.par
	for ((j=$smax1; j<=$smax2; j+=$smaxincr))
	do
		sed -E -i "s/[0-9]+[ \t]+smax/$j\tsmax/" LTS.par
		./LTS_polystyrene
		sed -i 's/D/E/g' LTS.sol
		cp LTS.sol LTS_ns${i}_smax$j.sol #necessarie le {}
	done
done

