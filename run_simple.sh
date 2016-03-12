#!/bin/bash
for i in $(seq 10)
do
	./ILS-UrApHMP 6 5 $i < instances/AP150.txt >> logs/AP150/AP150-6-5.csv
	echo "6 5: $i...ok!"
done
for p in 7 8
do
	for r in 2 3 4 5 6 7
	do
		if [ $r -eq $p ]
		then
			break
		fi
		for i in $(seq 10)
		do
			./ILS-UrApHMP $p $r $i < instances/AP150.txt >> logs/AP150/AP150-$p-$r.csv
			echo "$p $r: $i...ok!"
		done
	done
done
