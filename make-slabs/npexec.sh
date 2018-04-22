#!/bin/bash

if [ -e npslabs/ ]; then
	python3 nonpolar.py
else
	mkdir npslabs/
	python3 nonpolar.py
fi

cd npslabs/

for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 

do

i="d"$j

if [ -e $i ]; then
	a=$i/POSCAR.vasp
else
	mkdir $i
	a=$i/POSCAR.vasp
fi

printf "SYSTEM is ZnO Slabs\n1.0\n" > $a
cat $i"vecs" >> $a
printf "Zn O\n" >> $a
cat $i"num" >> $a
printf " " >> $a
cat $i"num" >> $a
printf "\nCartesian\n" >> $a
cat $i"zn1" >> $a
cat $i"zn2" >> $a
cat $i"o1" >> $a
cat $i"o2" >> $a

done

