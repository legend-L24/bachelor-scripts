#!/usr/bin/bash
./get_data2.py
let num_0=300
for i in $(seq 2 1 6)
do
	let num=num_0+100*i
        #mkdir ${i}
	cp 4.data ${i}/1.data
	cp lammps.sbs ${i}
	cp geom_${num}.in ${i}/geom.in
	echo ${num}
	#cp gap/* ${i}
	cd ${i}
	sbatch lammps.sbs
	cd ..
done
