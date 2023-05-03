#!/home/ytli/anaconda3/envs/quip/bin/python
from ase.io import read, write
import numpy as np
from ase.io.lammpsdata import read_lammps_data
atoms = [read_lammps_data("1.data", style="atomic")]
#atoms = read("1/geom.dump",":100:5000")
#atoms = read("mat.cif","::")
#atoms = read("1.data","::")
for atom in atoms:
    #if atom.symbol!='Na': continue
    print("check one structure")
    atom = atom[[i for i, j in enumerate(atom) if j.symbol!='Na']]
    dist_arr = atom.get_all_distances(mic=True)
    dist_arr += 5*np.eye(dist_arr.shape[0])
    if np.min(dist_arr)<1:
        print("there is still error",np.min(dist_arr))
