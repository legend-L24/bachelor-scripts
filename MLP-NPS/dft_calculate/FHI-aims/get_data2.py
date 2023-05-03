#!/u/ytli/anaconda3/envs/quip/bin/python
import os
import ase
from ase.io.lammpsdata import write_lammps_data
from ase.io import read, write
import numpy as np
import random as r

#config_file  = "train.xyz"
filelist = os.listdir(".")
atoms_kind = [ 'Na', 'P', 'S']
atoms = read("tz9b00322_si_001.cif", "::")
#atoms = read("tz9b00322_si_001.cif","::")
print(len(atoms[0]))
count = 1
for i in range(0, len(atoms)):
    #if i > 3: print("error!!"); continue
    atom = atoms[i]
    print("one time", i)
    cs = atoms[i].get_chemical_symbols()
    for idx, c in enumerate(cs):
        if c == atoms_kind[0]:
           cs[idx] = 'Ar'
        if c == atoms_kind[1]:
           cs[idx] = 'B'
        if c == atoms_kind[2]:
           cs[idx] = 'C'
    atom.set_chemical_symbols(cs)
    atom = atom.repeat([1,1,1])
    write_lammps_data(str(count)+'.data', atom, atom_style='atomic')
    count +=1
