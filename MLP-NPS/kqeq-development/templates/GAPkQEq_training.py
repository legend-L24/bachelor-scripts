import sys
sys.path.append("/datavon1/kQeqDevelopment/gradSOAP/kqeq/")

import random
import matplotlib.pyplot as plt

from kqeq.qeq import charge_eq
from kqeq.kernel import kernel
from kqeq.kernelGAPkQEq import soapGAPkQEq
from kqeq.kqeq import kernel_qeq
from kqeq.GAPkQEq import GAPkQEq#, get_energies_perAtom,get_energies
from kqeq.funct import  get_dipoles, get_charges, get_energies, get_energies_perAtom, create_lammps

import ase
from ase.units import Hartree
from ase import Atoms
from ase.io import read, write

import numpy as np
from sklearn.metrics import mean_squared_error
import math

# setting up energy of atoms in vaccum, and hardness (eta) 
atom_energy = {"Zn": -49117.02929728,
               "O": -2041.3604,
               "H": -13.63393,
               "Li": -203.1890559032085}


hard_lib = {"O": 1, "Zn": 1, "H": 0, "Li":0}

# load molecules, shuffle them and select certain number of structures
#all_molecules = read("/datavon1/DatabasesMolecules/ZnOData/ZnOall.xyz@:",format="extxyz")

#all_molecules = read("/datavon1/DatabasesMolecules/neutralMinima/ClusterData3_8.xyz@:",format="extxyz")
all_molecules = read("/datavon1/DatabasesMolecules/LiH/LiH.xyz@:",format="extxyz")
#all_molecules.extend(read("/datavon1/DatabasesMolecules/expandingWater/expand20/expHir.xyz@:",format="extxyz"))
#all_molecules.extend(read("/datavon1/DatabasesMolecules/expandingWater/expand16/expHir.xyz@:",format="extxyz"))
random.seed(20)
random.shuffle(all_molecules)

train_size = 100
val_size = 0
test_size = 50

training_set = all_molecules[:train_size]
valid_molecules = all_molecules[train_size:train_size+val_size]
testing_set = all_molecules[val_size+train_size:val_size+train_size+test_size]

# set up kernel
atsize = 1
rcut = 3.0
 


# Energy of the training set
ref_en_train = get_energies(mols=training_set,atom_energy = atom_energy)
desdictGAP = f"soap cutoff=5.0 zeta=2 delta = 1.0 atom_sigma =1.0 l_max=8 n_max=10 n_Z=2 Z={{1 3}} n_species=2 species_Z={{1 3}}"
desdictkQEq = f"soap cutoff=2.5 zeta=2 delta = 1.0 atom_sigma =0.5 l_max=8 n_max=10 n_Z=2 Z={{1 3}} n_species=2 species_Z={{1 3}}"

#desdictGAP = f"soap cutoff=5 zeta=2 delta = 1.0 atom_sigma =0.6 l_max=8 n_max=10 n_Z=2 Z={{1 8}} n_species=2 species_Z={{1 8}}"
#desdictkQEq = f"soap cutoff=2.5 zeta=2 delta = 1.0 atom_sigma =0.3 l_max=8 n_max=10 n_Z=2 Z={{1 8}} n_species=2 species_Z={{1 8}}"

SOAP_Kernel = soapGAPkQEq(multi_SOAP=False,
                        descriptor_dict_GAP=desdictGAP,
                        descriptor_dict_kQEq=desdictkQEq,
                        training_set=training_set)
# Create an instance of the kQEq class
my_kqeq = GAPkQEq(Kernel=SOAP_Kernel,
                     scale_atsize=atsize,
                     radius_type="qeq",
                     sparse=True,
                     sparse_count=1000,
                     hard_lib=hard_lib)

'''
my_kqeq.train(lambda_reg = 0.1, 
              lambda_charges = 10,
              lambda_charges_min= 0.0000001,
              atom_energy=atom_energy,
              iter_charges=30,
              charge_keyword='hirshfeld')
'''
my_kqeq.train(lambda_reg = 0.00001, 
            lambda_reg_grad = 1.2,
            lambda_charges = 0.01,
            lambda_charges_grad=3,
            atom_energy=atom_energy,
            iter_charges=5,
            charge_keyword='hirshfeld')

E_test = []
q_test = []
for a in testing_set:
    E,q,_ = my_kqeq.calculate(a)
    E_test.append(E/len(a))
    q_test.extend(q)
ref_en_test = get_energies_perAtom(mols=testing_set,atom_energy = atom_energy)


plt.title(f"Energy test")
plt.plot(ref_en_test,ref_en_test,color="black")
plt.scatter(ref_en_test, np.array(E_test),color="blue",alpha = 0.5)
plt.tight_layout()
plt.xlabel('Reference DFT Energy per atom / eV',fontsize=14)
plt.ylabel('kQEq Energy per atom / eV',fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()

'''
c_ref_train = get_charges(training_set,charge_keyword="hirshfeld")
c_ref_test = get_charges(testing_set,charge_keyword="hirshfeld")


plt.figure(figsize=(8, 5), dpi=100)
plt.scatter(c_ref_test,q_test,alpha = 0.5)
plt.plot(c_ref_test,c_ref_test,color = "black")
#plt.xlim([-8,8])
#plt.ylim([-8,8])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("test charger",fontsize=26)
plt.xlabel('Reference charge',fontsize=22)
plt.ylabel('kQEq charge',fontsize=22)
#plt.savefig("trainDip.pdf",bbox_inches='tight')
plt.show()


E_train = []
q_train = []
for a in training_set:
    E,q,_ = my_kqeq.calculate(a)
    E_train.append(E/len(a))
    q_train.extend(q)
ref_en_train = get_energies_perAtom(mols=training_set,atom_energy = atom_energy)

plt.title(f"Energy train")
plt.plot(ref_en_train,ref_en_train,color="black")
plt.scatter(ref_en_train, E_train,color="blue",alpha = 0.5)
plt.tight_layout()
plt.xlabel('Reference DFT Energy per atom / eV',fontsize=14)
plt.ylabel('kQEq Energy per atom / eV',fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='minor', labelsize=10)
plt.show()

plt.figure(figsize=(8, 5), dpi=100)
plt.scatter(c_ref_train,q_train,alpha = 0.5)
plt.plot(c_ref_train,c_ref_train,color = "black")
#plt.xlim([-8,8])
#plt.ylim([-8,8])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("train charger",fontsize=26)
plt.xlabel('Reference charge',fontsize=22)
plt.ylabel('kQEq charge',fontsize=22)
#plt.savefig("trainDip.pdf",bbox_inches='tight')
plt.show()
'''
print("DOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOONE")
my_kqeq.save_model(atom_energy=atom_energy)
E = 0
for atom in testing_set[0]:
    E += atom_energy[atom.symbol]
_,q,eneg = my_kqeq.calculate(testing_set[0])
E_all,E_gap,E_kqeq = my_kqeq.calculateEnergy(testing_set[0])
print("Without nonatom E: ", E)
print("GAP: ", E_gap)
print("kQeq: ", E_kqeq)
print(E+E_gap+E_kqeq)
print(E_all+E)
print("ref: ",testing_set[0].info["energy"])
write("test.xyz",testing_set[0])
create_lammps(testing_set[0],q,eneg,["H","Li"],None,"hello",hard_lib)

write("testing_set.xyz",testing_set)
write("training_set.xyz",training_set)