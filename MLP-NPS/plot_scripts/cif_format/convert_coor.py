#!/home/ytli/anaconda3/envs/quip/bin/python

from ase.io import read
import sys
import json
from numpy.linalg import norm
from numpy import argmin
json_path = "/home/ytli/decryst_cp/python/poly.json"
ligand = "L7"

filename = './'+sys.argv[1]
poly = read(filename)
del poly[[atom.index for atom in poly if atom.symbol=='H']]

poly.translate(-poly.get_center_of_mass())

poss_norm = norm(poly.get_positions(),axis=1)
center_atom = poly.pop(argmin(poss_norm))
print("this is center_atom: ", center_atom.symbol)    
poly.translate(-center_atom.position)

# maybe you can remove H atom to reduce problems
new_data = {ligand: {}}
new_data[ligand]["elems"] = dict(map(reversed, enumerate(set(poly.get_chemical_symbols()))))
new_data[ligand]["region"] = "0 180 0 180 0 180"
new_data[ligand]["ref_coor"] = [[i] + j for i, j in zip(poly.get_chemical_symbols(),[list(i) for i in poly.get_positions()])]
with open(json_path, "r") as f:
    file = f.read()
    if len(file) > 0:
        old_data = json.load(open(json_path))
    else: 
        old_data = {}
    old_data.update(new_data)

with open(json_path, "w") as f:
    data = json.dumps(old_data, indent=1)
    f.write(data)

