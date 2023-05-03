#!/u/ytli/anaconda3/envs/quip/bin/python
import numpy as np
from ase.io import read
import os

path_ls = ["geom.dump"]#, "./out_2.lammpstrj"]
pick_ls = ["-100::10"]# "17:38:3"]
basis_path = "/nexus/posix0/FHI-Theory/Software/FHI-aims/species_defaults/defaults_2020/intermediate/"
num=0
atoms_kind = ['Na', 'P', 'S']
for path, pick in zip(path_ls, pick_ls):
    atoms = read(path,"::")
    print("this is the whole legenth:", len(atoms))
    atoms = read(path, pick)
    for atom in atoms:
        num += 1
        cs = atom.get_chemical_symbols()
        for jdx, c in enumerate(cs):
            if c == "H":
                cs[jdx] = atoms_kind[0]
            if c == "He":
                cs[jdx] = atoms_kind[1]
            if c == "Li":
                cs[jdx] = atoms_kind[2]
        atom.set_chemical_symbols(cs)
        cell = atom.get_cell()
        pos = atom.get_positions()
        w_str = "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[0][0], cell[0][1], cell[0][2])
        w_str += "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[1][0], cell[1][1], cell[1][2])
        w_str += "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[2][0], cell[2][1], cell[2][2])

        w_str += "\n\n"
        for idx, ps in enumerate(pos):
            w_str += "atom %12.6f %12.6f %12.6f %2s\n" %(ps[0], ps[1], ps[2], cs[idx])
        with open("geometry_"+str(num)+".in", "w") as f:
            f.write(w_str)
		#print(w_str)
		#this is control.in
        w_str = """
		    xc              pbe
		    relativistic    atomic_zora scalar
		    k_grid          4 4 4
			sc_accuracy_forces  0.0001
		"""
        for file in os.listdir(basis_path):
            for elem in set(cs):
                if f"_{elem}_" in file:
                    with open(basis_path+file) as f:
                        w_str += f.read()
        with open("control_"+str(num)+".in", "w") as f:
            f.write(w_str)
print(num)
