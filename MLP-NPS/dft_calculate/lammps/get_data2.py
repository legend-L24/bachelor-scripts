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
config_file = "../na3ps4_0315.xyz"
atoms_kind = [ 'Na', 'P', 'S']
atoms = read(config_file, "-10::")
cs = atoms[0].get_chemical_symbols()
for idx, c in enumerate(cs):
    if c == atoms_kind[0]:
       cs[idx] = 'Ar'
    if c == atoms_kind[1]:
       cs[idx] = 'B'
    if c == atoms_kind[2]:
       cs[idx] = 'C'

count = 1
for i in np.arange(1, len(atoms)):
    #if i > 3: print("error!!"); continue
    atom = atoms[i]
    print(len(atom))
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
