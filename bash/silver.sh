#!/bin/bash

#funzzziona!
pos1=$1
pos3=$2

./txt_to_dat.R $pos1 $pos3

#errore
in=$(echo $pos1 | sed 's/_pos1\.txt$/.dat/')
cp EXT_$in LTS.dat
./LTS_silver
sed -i 's/D/E/g' LTS.sol > DIST_$in

out=$(echo $pos1 | sed 's/pos1\.txt$/.png/')
gnuplot -e "filename = '${out}'" dist.gp

#errlog
#cp ext2.dat LTS.dat
#./LTS_polystyrene
#sed -i 's/D/E/g' LTS.sol

#out2=$(echo $pos1 | sed 's/pos1\.txt$/\errlog.png/')
#gnuplot -persist -e "filename = '${out2}'" dist.gp




