import sys
sys.path.append("/datavon1/kQeqDevelopment/kQeqJax/kqeq/")


from matplotlib.ticker import LinearLocator
import time
import random
import matplotlib.pyplot as plt

from kqeq.qeq import charge_eq
from kqeq.kernel import SOAPKernel
from kqeq.kQEq import kernel_qeq
from kqeq.funct import *

import ase
from ase.units import Hartree
from ase import Atoms
from ase.io import read, write

import numpy as np
from sklearn.metrics import mean_squared_error
import math

dataset = read("/datavon1/DatabasesMolecules/neutralMinima/ClusterData3_8.xyz@:",format="extxyz")

train_set, val_set, test_set = prep_structures(dataset,100,20, 50)

atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*Hartree}

hard_lib = {"O": 1, "Zn": 1, "H": 0, "Li":0}

desdict = {"nmax": 4,
           "lmax": 3,
           "rcut": 3.0,
           "sigma": 3.0 / 6,
           "periodic": False}

SOAP_Kernel = SOAPKernel(multi_SOAP=False,
                     descriptor_dict=desdict,
                     training_set=train_set,
                     training_system_charges=[0 for a in train_set],
                     validation_set=None)

my_kqeq = kernel_qeq(Kernel=SOAP_Kernel,
                     scale_atsize=1.0,
                     radius_type="qeq",
                     sparse=True,
                     sparse_count=1000,
                     hard_lib=hard_lib)

my_kqeq.train(target_lambdas =[0.01],#,0.01], 
              targets = ["dipole"],#,"charges"], 
              charge_keyword="hirshfeld")

kqeq_dipoles_train,kqeq_charges_train, _ = my_kqeq.predict(train_set)
kqeq_dipoles_test,kqeq_charges_test, _ = my_kqeq.predict(test_set)
d_ref_train = get_dipoles(train_set)#*2.541746229
d_ref_test = get_dipoles(test_set)#*2.541746229
c_ref_train = get_charges(train_set,charge_keyword="hirshfeld")
c_ref_test = get_charges(test_set,charge_keyword="hirshfeld")

plot_basics(kqeq_dipoles_test,d_ref_test,preset="dipole",save="hokus")
plot_adv(kqeq_charges_test,c_ref_test,preset="charges")
