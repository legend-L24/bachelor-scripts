#!/home/ytli/anaconda3/envs/quip/bin/python
from ase.io import read, write
import numpy as np
import sys
#atoms = read("train.xyz","::")
atoms = read(sys.argv[1], index="::", format="cif")
for atom in atoms:
    print("check one structure")
    dist_arr = atom.get_all_distances(mic=True)
    dist_arr += 5*np.eye(dist_arr.shape[0])
    if np.min(dist_arr)<1.5:
        print("there is still error", np.min(dist_arr))
