#!/u/ytli/anaconda3/envs/lyt/bin/python
import numpy as np
from ase.io import read, write
import os

basis_path = "/nexus/posix0/FHI-Theory/Software/FHI-aims/species_defaults/defaults_2020/light/" 
atoms = read("/u/ytli/sodium_train.xyz", "::")
for jdx, atom in enumerate(atoms):
    cell = atom.get_cell()
    pos = atom.get_positions()
    cs = atom.get_chemical_symbols()
    w_str = "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[0][0], cell[0][1], cell[0][2])
    w_str += "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[1][0], cell[1][1], cell[1][2])
    w_str += "lattice_vector %12.6f %12.6f %12.6f \n" %(cell[2][0], cell[2][1], cell[2][2])

    w_str += "\n\n"
    for idx, ps in enumerate(pos):
        w_str += "atom %12.6f %12.6f %12.6f %2s\n" %(ps[0], ps[1], ps[2], cs[idx])
    with open("geometry_"+str(jdx)+".in", "w") as f:
        f.write(w_str)
    #print(w_str)
    
    #this is control.in
    w_str = """
        xc              pbe
        relativistic    atomic_zora scalar
        k_grid          1 1 1
        relax_geometry  bfgs 1e-1
    """
    for file in os.listdir(basis_path):
        for elem in set(cs):
            if f"_{elem}_" in file:
                with open(basis_path+file) as f:
                    w_str += f.read()
    with open("control_"+str(jdx)+".in", "w") as f:
        f.write(w_str)
