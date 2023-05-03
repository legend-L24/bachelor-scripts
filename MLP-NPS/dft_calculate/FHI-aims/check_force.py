#!/u/ytli/anaconda3/envs/quip/bin/python

from ase.io import read, write
import numpy as np
atoms = read("gama.xyz", "::")
atoms_new = []
print("this is original number", len(atoms))
for idx, atom in enumerate(atoms):
    frc_arr = atom.get_array('forces')
    print(np.max(frc_arr),np.min(frc_arr))
    if np.max(frc_arr)<4 and np.min(frc_arr)>-4:
        atoms_new.append(atom)
    else: print("this is the wrong one",idx)
print("this is valid number: ", len(atoms_new))
write("result.xyz", atoms_new)
'''
frc_gap = np.array([])
ener_gap = []
for idx, atom in enumerate(atoms):
    frc_gap = np.append(frc_gap, atom.get_array('force'))
    atom_len = len(atom.get_chemical_symbols())
    ener_dft.append(atoms_[idx].info['energy']/atom_len)
    ener_gap.append(atom.get_potential_energy()/atom_len)
ener_gap = np.array(ener_gap)
for i in frc_gap:
'''
