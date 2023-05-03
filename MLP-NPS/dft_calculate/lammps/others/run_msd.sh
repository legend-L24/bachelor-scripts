#!/usr/bin/bash
./get_data2.py
for i in $(seq 1 1 5)
do
	#mkdir ${i}
	cp ../read_msd.py ${i}
	cd ${i}
	./read_msd.py
	cd ..
done
