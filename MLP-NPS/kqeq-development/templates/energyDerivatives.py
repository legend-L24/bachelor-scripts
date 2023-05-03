from ase.io import read,write
import sys
sys.path.append("/datavon1/kQeqDevelopment/GAPforces/kqeq")

import numpy as np
import random

from kqeq.kernel import SOAPKernel
from kqeq.funct import *
import matplotlib.pyplot as plt
from ase.units import Bohr, Hartree

atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*Hartree}

def num_grad(mol,my_kqeq,h=0.0001,direction=0,iatom=0):
    tmpmol = mol.copy()
    pos = tmpmol.get_positions()#/Bohr
    pos[iatom][direction] += h
    tmpmol.set_positions(pos)#*Bohr)
    res = my_kqeq.calculate(tmpmol)
    Eplus = res['energy']#/Hartree
    pos[iatom][direction] += -2.0*h
    tmpmol.set_positions(pos)#*Bohr)
    res = my_kqeq.calculate(tmpmol)
    Eminus = res['energy']#/Hartree
    pos[iatom][direction] += h
    tmpmol.set_positions(pos)#*Bohr)
    energy = (Eplus-Eminus)/(2.0*h)#*Bohr)
    return energy#*Hartree#/Bohr
    
def num_grads(mol,my_kqeq,h=0.0001):
    f = np.zeros((len(mol),3))
    for i in range(len(mol)):
        for direction in range(3):
            f[i,direction] = -num_grad(mol,my_kqeq,h=h,direction=direction,iatom=i)
    return f

all_molecules = read("/datavon1/DatabasesMolecules/LiH/LiH.xyz@:",format="extxyz")
random.seed(20)
random.shuffle(all_molecules)
val_size = 10
train_size = 100
test_size = 50
train_set = all_molecules[:train_size]
valid_molecules = all_molecules[train_size:train_size+val_size]
test_set = all_molecules[val_size+train_size:val_size+train_size+test_size]

from kqeq.GAP import GAP
desdictGAP = {"nmax": 6,"lmax": 5,"rcut": 5.0,"sigma": 5.0/8,"periodic": False}
SOAP_Kernel = SOAPKernel(multi_SOAP=False,descriptor_dict=desdictGAP,training_set=train_set,training_system_charges=[0 for a in train_set],validation_set=None)
ref_en_train = get_energies(mols=train_set,atom_energy = atom_energy)
my_kqeq = GAP(Kernel=SOAP_Kernel, sparse=True, sparse_method = "CUR", sparse_count=1000)
my_kqeq.train(atom_energy=atom_energy, lambda_reg = 0.01)


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
f_num = num_grads(test_set[2],my_kqeq)
results = my_kqeq.calculate(test_set[2])
f = results['forces']
print(f)
print("$$$$$$$$$$$$$$$$$$")
print(f_num)
print("$$$$$$$$$$$$$$$$$$")
print(f_num/f)
''' KQEq

atsize, radtype, rcut = 1.0, "qeq", 3.5
from kqeq.kqeq import kernel_qeq
from kqeq.kernelgap import kernelGAP
desdictkQEq = f"soap cutoff=2.5 zeta=2 delta = 1.0 atom_sigma =0.5 l_max=8 n_max=10 n_Z=2 Z={{1 3}} n_species=2 species_Z={{1 3}}"
SOAP_Kernel = kernelGAP(multi_SOAP=False,descriptor_dict=desdictkQEq,training_set=training_set,validation_set=valid_molecules)
#from kqeq.kernel import kernel
desdict = {"nmax": 6, "lmax": 3, "rcut": rcut, "sigma": rcut / 8, "periodic": False}
#SOAP_Kernel = kernel(Kernel='SOAP',Descriptor='SOAP',multi_SOAP=False,descriptor_dict=desdict,training_set=training_set,validation_set=valid_molecules)
my_kqeq = kernel_qeq(Kernel=SOAP_Kernel,scale_atsize=atsize,radius_type=radtype,sparse=True,sparse_method="CUR",calculate_kernel_forces=True)
my_kqeq.train(target_lambdas=[0.1], targets = ["dipole"])
#print(my_kqeq.weights)
print(len(testing_set))
f_num = num_grads(testing_set[2],my_kqeq)
results = my_kqeq.calculate(testing_set[2])
f = results['forces']

print(f)
print("$$$$$$$$$$$$$$$$$$")
print(f_num)
print("$$$$$$$$$$$$$$$$$$")
print(f_num/f)

desdict = {"nmax": 6, "lmax": 3, "rcut": rcut, "sigma": rcut / 8, "periodic": False}
SOAP_Kernel = kernel(Kernel='SOAP',Descriptor='SOAP',multi_SOAP=False,descriptor_dict=desdict,training_set=training_molecules,validation_set=valid_molecules)

my_kqeq.train(lambda_reg=0.1, wq = 0.5)

results = my_kqeq.calculate(testing_set[1])
f = results['forces']

f_num = num_grads(testing_set[1],my_kqeq)
print(f)
print(f_num)
#print(f/f_num)
#print(f_num/f)
'''