#!/bin/bash

if [ -e pslabs/ ]; then
	python3 polar.py
else
	mkdir pslabs/
	python3 polar.py
fi

cd pslabs/

for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19

do

i="d"$j
if [ -e $i ]; then 
	a=$i/POSCAR.vasp
else
	mkdir $i
	a=$i/POSCAR.vasp
fi

if [ -e ads_id ]; then
	printf "SYSTEM is ZnO Slabs with adsorbate\n1.0\n" > $a
	cat $i"vecs" >> $a
	printf "Zn O " >> $a
	cat ads_id | tr '\n' ' ' >> $a
	printf "\n" >> $a
	cat $i"num" >> $a
	printf " " >> $a
	cat $i"num" >> $a
	printf " " >> $a
	cat $i"tadsnum" >> $a
	printf " " >> $a
	cat $i"badsnum" >> $a
	printf "\nCartesian\n" >> $a
	cat $i"zn1" >> $a
	cat $i"zn2" >> $a
	cat $i"o1" >> $a
	cat $i"o2" >> $a
	cat $i"tads" >> $a
	cat $i"bads" >> $a
else
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

fi

done

