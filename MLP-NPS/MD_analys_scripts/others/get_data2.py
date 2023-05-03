#!/u/ytli/anaconda3/envs/quip/bin/python
import os
import ase
from ase.io.lammpsdata import write_lammps_data
from ase.io import read, write
import numpy as np
import random as r

#config_file  = "train.xyz"
filelist = os.listdir(".")
ciflist = []
for file in filelist:
    if file[-3:] == "cif": ciflist.append(file)
print(ciflist)
config_file = ciflist[1]
config_file = "Na4P2S6.xyz"
atoms_kind = [ 'Na', 'P', 'S']
atoms = read(config_file, "::")
cs = atoms[0].get_chemical_symbols()
for idx, c in enumerate(cs):
    if c == atoms_kind[0]:
       cs[idx] = 'Ar'
    if c == atoms_kind[1]:
       cs[idx] = 'B'
    if c == atoms_kind[2]:
       cs[idx] = 'C'
atoms = read("tz9b00322_si_001.cif","::")
print(len(atoms[0]))
count = 1
for i in range(0, len(atoms)):
    #if i > 3: print("error!!"); continue
    atom = atoms[i]
    cs = atoms[i].get_chemical_symbols()
    for idx, c in enumerate(cs):
        if c == atoms_kind[0]:
           cs[idx] = 'Ar'
        if c == atoms_kind[1]:
           cs[idx] = 'B'
        if c == atoms_kind[2]:
           cs[idx] = 'C'
    atom.set_chemical_symbols(cs)
    atom = atom.repeat([2,1,1])
    write_lammps_data(str(count)+'.data', atom, atom_style='atomic')
    count +=1
