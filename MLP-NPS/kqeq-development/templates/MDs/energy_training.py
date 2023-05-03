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

#dataset = read("/datavon1/DatabasesMolecules/neutralMinima/ClusterData3_8.xyz@:",format="extxyz")
#dataset = read("/home/yli/dataset/LiH.xyz@:",format="extxyz")
'''
dataset = read("/home/yli/software/kqeq-development/templates/MDs/zero.xyz@:",format="extxyz")
train_set, val_set, test_set = prep_structures(dataset,100,20, 50)
'''

dataset = read("/home/yli/ZnO_filt.xyz@:",format="extxyz")
#dataset = read("/home/yli/dataset/LiH.xyz@:",format="extxyz")
train_set, val_set, test_set = prep_structures(dataset,41,0, 0)

'''
dataset = read("/home/yli/software/kqeq-development/templates/MDs/zero.xyz@:",format="extxyz")
train_set, val_set, test_set = prep_structures(dataset,100,20, 50)
'''
atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*Hartree}


hard_lib = {"O": 1, "Zn": 1, "H": 0, "Li":0}



desdictSOAP = {"nmax": 6,
           "lmax": 5,
           "rcut": 3.0,
           "sigma": 3.0/8,
           "periodic": False}


SOAP_Kernel = SOAPKernel(multi_SOAP=False,
                     descriptor_dict=desdictSOAP,
                     training_set=train_set,
                     training_system_charges=[0 for a in train_set],
                     validation_set=None)

# Energy of the training set
ref_en_train = get_energies(mols=train_set,atom_energy = atom_energy)

# Create an instance of the kQEq class
my_kqeq = kernel_qeq(Kernel=SOAP_Kernel,
                     scale_atsize=1.0,
                     radius_type="qeq",
                     sparse=True,
                     sparse_count=1000,
                     hard_lib=hard_lib)

my_kqeq.trainEnCharge(atom_energy=atom_energy,
                      charge_keyword="hirshfeld", 
                      lambda_reg = 0.01,
                      lambda_target = 0.2,
                      lambda_target_min=0.001,
                      iter_charges = 5)

E_train = []
for a in train_set:
    E_train.append(my_kqeq.calculate(a)["energy"]/len(a))
ref_en_train = get_energies_perAtom(mols=train_set,atom_energy = atom_energy)

E_test = []
for a in test_set:
    E_test.append(my_kqeq.calculate(a)["energy"]/len(a))
ref_en_test = get_energies_perAtom(mols=test_set,atom_energy = atom_energy) 

plot_basics(E_train,ref_en_train,preset="energy", save="energy_train_1.png")
#plot_basics(E_test,ref_en_test,preset="energy", save="energy_test_1.png")

#my_kqeq.save_kQEq()
#write("md_beg.xyz",test_set[10])
#write('training_set.xyz',train_set)


