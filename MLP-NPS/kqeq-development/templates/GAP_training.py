import sys
sys.path.append("/datavon1/kQeqDevelopment/GAPforces/kqeq/")

import random
import matplotlib.pyplot as plt

from kqeq.GAP import GAP
from kqeq.funct import *
from kqeq.kernel import SOAPKernel

import ase
from ase.units import Hartree
from ase import Atoms
from ase.io import read, write

import numpy as np
from sklearn.metrics import mean_squared_error
import math

#dataset = read("/datavon1/DatabasesMolecules/neutralMinima/ClusterData3_8.xyz@:",format="extxyz")
dataset = read("/datavon1/DatabasesMolecules/LiH/LiH.xyz@:",format="extxyz")
train_set, val_set, test_set = prep_structures(dataset,100,20, 50)

atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*Hartree}

hard_lib = {"O": 1, "Zn": 1, "H": 0, "Li":0}



desdictGAP = {"nmax": 6,
           "lmax": 5,
           "rcut": 5.0,
           "sigma": 5.0/8,
           "periodic": False}


SOAP_Kernel = SOAPKernel(multi_SOAP=False,
                     descriptor_dict=desdictGAP,
                     training_set=train_set,
                     training_system_charges=[0 for a in train_set],
                     validation_set=None)
# Energy of the training set
ref_en_train = get_energies(mols=train_set,atom_energy = atom_energy)

# Create an instance of the kQEq class
my_kqeq = GAP(Kernel=SOAP_Kernel,
                     sparse=True,
                     sparse_method = "CUR",
                     sparse_count=1000)

my_kqeq.train(atom_energy=atom_energy, 
                      lambda_reg = 0.01)



E_train = []
for a in train_set:
    E_train.append(my_kqeq.calculate(a)["energy"]/len(a))
ref_en_train = get_energies_perAtom(mols=train_set,atom_energy = atom_energy)

E_test = []
for a in test_set:
    E_test.append(my_kqeq.calculate(a)["energy"]/len(a))
ref_en_test = get_energies_perAtom(mols=test_set,atom_energy = atom_energy) 

plot_basics(E_train,ref_en_train,preset="energy")
plot_basics(E_test,ref_en_test,preset="energy")





