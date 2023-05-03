import sys
sys.path.append("/home/yli/software/kqeq-development")

import random
import matplotlib.pyplot as plt

from kqeq.qeq import charge_eq
from kqeq.kernel import *
from kqeq.kQEq import kernel_qeq
from kqeq.funct import *

import ase
from ase.units import Hartree
from ase import Atoms
from ase.io import read, write

import numpy as np
from sklearn.metrics import mean_squared_error
import math
import matplotlib.pyplot as plt
#dataset = read("/datavon1/DatabasesMolecules/neutralMinima/ClusterData3_8.xyz@:",format="extxyz")
#dataset = read("/home/yli/dataset/LiH.xyz@:",format="extxyz")
'''
from ase.io import read
atoms = read("/home/yli/ZnO_0.xyz","::")
pos = []
for atom in atoms:
    ls = list(atom.arrays["hirshfeld"])
    ls = [float(i) for i in ls]
    atom.set_array("hirshfeld", np.array(ls, dtype=float))
    pos.append(atom)
print(type(pos[0].arrays["hirshfeld"][0]))
write("/home/yli/ZnO.xyz", pos)

dataset = read("/home/yli/ZnO.xyz@:",format="extxyz")
q_arr = dataset[0].arrays["hirshfeld"]
print(type(q_arr[0]))
'''
dataset = read("/home/yli/ZnO_filt.xyz@:",format="extxyz")
dataset = read("/home/yli/software/kqeq-development/templates/MDs/zero.xyz@:",format="extxyz")
q_arr = dataset[-2].arrays["hirshfeld"]
#print(q_arr)
q_arr = q_arr.astype(float)
#print(type(q_arr))
qeq_1 = charge_eq(dataset[-2],q=q_arr)
#qeq_1.calc_Eqeq_old()
print(qeq_1.E) 
qeq_1.calc_Eqeq()
print(qeq_1.E)

#H = qeq_1.get_H()
#Abar = qeq_1.get_Abar()
#np.linalg.inv(H)
#print(H)
'''
atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*Hartree}
ref_ener = atom_energy["Zn"]*108 + atom_energy["O"]*108
ener_ls = [(i.get_potential_energy()-ref_ener)/216 for i in dataset]
plt.plot(ener_ls)
print(ener_ls)
idx_ls = []
for idx, ener in enumerate(ener_ls):
    if ener<-3.8: idx_ls.append(idx)
print(len(idx_ls), len(ener_ls))

pos = []
for idx in idx_ls:
    pos.append(dataset[idx])
write("/home/yli/ZnO_filt.xyz", pos)
plt.savefig("/home/yli/ener.png")
'''
