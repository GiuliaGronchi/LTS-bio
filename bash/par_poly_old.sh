#!/bin/bash

nsmin=$1
nsmax=$2
incr=$3
smax1=$4
smax2=$5
sincr=$6

for ((j=$smax1; j<=$smax2; j+=$sincr))
do
	sed -E -i "s/[0-9]+[ \t]+smax/$j\tsmax/" LTS.par
	for ((i=$nsmin; i<=$nsmax; i+=$incr))
		do
		sed -E -i "s/[0-9]+[ \t]+ns/$i\tns/" LTS.par
		./LTS_polystyrene
		sed -i 's/D/E/g' LTS.sol
		cp LTS.sol LTS_ns${i}_smax$j.sol
	done
done
