#!/usr/bin/bash

for i in $(seq 3 1 10)
do
	cd ${i}
	sbatch gap.sbs
	cd ..
done
