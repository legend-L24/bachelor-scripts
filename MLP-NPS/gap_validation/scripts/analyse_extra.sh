#!/usr/bin/bash

#cp val_1.xyz test.xyz
#/home/ytli/softwares/QUIP/build/linux_x86_64_gfortran_openmp/quip E=T F=T atoms_filename=test.xyz param_filename=gp_2b_soap.xml | grep AT | sed 's/AT//' > quip_train.xyz
#../error.py

cp val_1.xyz test.xyz
/home/ytli/softwares/QUIP/build/linux_x86_64_gfortran_openmp/quip E=T F=T atoms_filename=test.xyz param_filename=gp_2b_soap.xml | grep AT | sed 's/AT//' > quip_train.xyz
../error_2.py
