#!/usr/bin/bash

for i in $(seq 3 1 10)
do
	cd ${i}
	nohup ../analyse_hpc.sh &
	cd ..
done
