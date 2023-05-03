import sys
sys.path.append("/datavon1/kQeqDevelopment/kQeqJax/kqeq/")
import time
import random
import matplotlib.pyplot as plt
from kqeq.qeq import charge_eq
from kqeq.kernel import SOAPKernel
from kqeq.kQEq import kernel_qeq
from kqeq.funct import get_whole_energies_perAtom, get_dipoles, get_charges, get_energies, get_energies_perAtom, save_charges, get_whole_energies
import ase
from ase import Atoms
from ase.io import read, write
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

train_set = read("training_set.xyz@:",format="extxyz")

atom_energy = {"Zn": -49117.02929728, "O": -2041.3604,
               "H": -13.63393, "Li": -7.467060138769*units.Hartree}

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
                     sparse_count=1000)
my_kqeq.load_kQEq()
main = read("md_beg.xyz")
from kqeq.calculator import Kqeq
calc = Kqeq(my_kqeq)


main.calc = calc



MaxwellBoltzmannDistribution(main, temperature_K=300)
logging_args = dict(trajectory='md_prdEn.traj', logfile='md_prdEn.log', loginterval=1)
dyn = VelocityVerlet(main, 5 * units.fs,**logging_args)  # 5 fs time step.

import time
t0 = time.time()
def printenergy(a):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    char = a.calc.results["charges"]
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
    #print(char)

# Now run the dynamics
printenergy(main)
for i in range(1):
    dyn.run(50)
    printenergy(main)
print(time.time()-t0)
