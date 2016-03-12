#!/bin/bash
for i in $(seq 30);
do
	./ILS-UrApHMP 5 3 $i < AP60.txt >> logs/ILS.csv
	echo "$i...ok!"
done;
