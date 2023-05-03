#!/usr/bin/bash
for i in $(seq 1 1 3)
do
	cd glass_tetr_${i}
	for j in $(seq 4 1 8)
	do
		cd ${j}00K
		cp ../1.data .
		sbatch < lammps.sbs
		echo $i and $j
		cd ..
	done
	cd ..
done
