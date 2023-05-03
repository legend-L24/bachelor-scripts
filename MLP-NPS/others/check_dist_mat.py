#!/home/ytli/anaconda3/envs/quip/bin/python
from ase.io import read, write
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter, CifParser
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.data.LammpsData import from_structure
#atoms = read("1/geom.dump",":100:5000")
atoms = read("cm1c01113_si_006.cif","::")
parser = CifParser("cm1c01113_si_006.cif")
structure = parser.get_structures()[0]
#from_structure(structure)
#lammps_data = CombinedData([structure])
#lammps_data.write_file("1.data")
#cif = CifWriter(structure)
#cif.write_file("1.data")
#cif.write_file("mat.cif")
'''
atoms = [AseAtomsAdaptor.get_atoms(structure)]
for atom in atoms:
    #if atom.symbol!='Na': continue
    print("check one structure")
    dist_arr = atom.get_all_distances(mic=True)
    dist_arr += 5*np.eye(dist_arr.shape[0])
    if np.min(dist_arr)<1:
        print("there is still error",np.min(dist_arr))
'''
