import sys
sys.path.append("/datavon1/kQeqDevelopment/kQeqJax/kqeq/")
import time
import random
import matplotlib.pyplot as plt
from kqeq.qeq import charge_eq
from kqeq.kernel import SOAPKernel
from kqeq.kQEq import kernel_qeq
from kqeq.forMD import kQEqMD
from kqeq.funct import * 
import ase
from ase import Atoms
from ase.io import read, write
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units


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
                     validation_set=None)


# Create an instance of the kQEq class
my_kqeq = kQEqMD(Kernel=SOAP_Kernel,
                     scale_atsize=1.0,
                     radius_type="qeq",
                     sparse=True)
my_kqeq.load_kQEq()

train_set = read("training_set.xyz@:",format="extxyz")

test_set = read("testing_set.xyz@:", format="extxyz")


main = read("md_beg.xyz")
from kqeq.calculator import Kqeq
calc = Kqeq(my_kqeq)

main.calc = calc

print(my_kqeq.Kernel.training_set)

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
for i in range(5):
    dyn.run(50)
    printenergy(main)
print(time.time()-t0)
