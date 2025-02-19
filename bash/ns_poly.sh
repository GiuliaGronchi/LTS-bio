#!/bin/bash

#Permette di variare automaticamente ns (numero di punti) in LTS.par

e1=$1
e2=$2
incr=$3

for ((i=$e1; i<=$e2; i+=$incr))
do
	sed -E -i "s/[0-9]+[ \t]+ns/$i\tns/" LTS.par
	./LTS_polystyrene
	sed -i 's/D/E/g' LTS.sol
	cp LTS.sol LTS$i.sol
done
