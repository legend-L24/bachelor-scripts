#!/home/ytli/anaconda3/envs/quip/bin/python
from ase.io import read, write
import numpy as np
#atoms = read("1/geom.dump",":100:5000")
atoms = read("mat.cif","::")
atoms = read("cm1c01113_si_006.cif","::")
print(len(atoms[0]))
from ase.io.lammpsdata import write_lammps_data
write_lammps_data('1.data', atoms[0], atom_style='atomic')
for atom in atoms:
    print("atomic number:",len(atom))
    #if atom.symbol!='Na': continue
    print("check one structure")
    #atom = atom[[i for i, j in enumerate(atom) if j.symbol!='Na']]
    dist_arr = atom.get_all_distances(mic=True)
    print(dist_arr.shape)
    dist_arr += 5*np.eye(dist_arr.shape[0])
    if np.min(dist_arr)<1:
        print("there is still error",np.min(dist_arr))
    for i in range(dist_arr.shape[0]):
        for j in range(dist_arr.shape[1]):
            if dist_arr[i,j]<0.20:
                print(atom[i].symbol,atom[j].symbol,dist_arr[i,j])
