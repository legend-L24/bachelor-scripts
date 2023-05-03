import random
from turtle import color
import numpy as np
from ase.io import *
from ase.data import covalent_radii
from ase.units import Bohr,Hartree
from kqeq.qeq import charge_eq
from kqeq.data import uff_radius_qeq

def prep_structures(dataset,NTrain, NVal, NTest,file_extension = "xyz", file_format = "extxyz",seed = 20):
    nmols = len(dataset)
    random.seed(seed) 
    random.shuffle(dataset)
    train_set = dataset[:NTrain]
    valid_set = dataset[NTrain:NTrain+NVal]
    test_set = dataset[NVal+NTrain:NVal+NTrain+NTest]
    
    return train_set, valid_set, test_set

def get_charges(mols,charge_keyword):
    """Helper function to get 1D array of atomic charges from list of atoms objects.

    Parameters
    ----------
        mols : list
            List of ASE atoms objects with initial charges
    Returns
    -------
        np.array
            Numpy array containing the charges
    """
    ref_qs = []
    for mol in mols:
        ref_qs.extend(mol.arrays[charge_keyword])    
    return np.array(ref_qs)
    
def save_charges(mols, charges,charge_keyword = "kqeq_charge",name = "test.xyz"):
    
    count = 0

    for mol in mols:
        ref_qs = charges[count:count+len(mol)]
        count += len(mol)
        mol.arrays[charge_keyword] = ref_qs     
    write(name,mols)

def get_dipoles(mols):
    """Helper function to get 1D array of dipole vector elements from list of atoms objects.

    Parameters
    ----------
        mols : list
            List of ASE atoms objects with 'dipole_vector' in info (atomic units!)
    Returns
    -------
        np.array
            Numpy array containing the dipole vector elements
    """
    ref_mus = []
    for mol in mols:
        ref_mus.extend(mol.info['dipole_vector'])
    
    return np.array(ref_mus)

def get_energies(mols, atom_energy, energy_keyword = "energy"):
    """Function to get atomization energy. This function is used for training on energy, so energy is 
    returned in Hartrees. Dictionary of energies of atoms in vacuum is needed.    

    Parameters
    ----------
        mols : list
            list of ASE atoms with energies
        atom_energy : dictionary
            dictionary of energies of atoms in vacuum
        energy_keyword : string
            word used for energy extraction (energies has to be in eV) 

    Returns
    -------
        np.array
            Numpy array containing the energies in Hartree
    """
    energy = []
    for mol in mols:
        mol_energy = mol.info[energy_keyword]
        for element in atom_energy:
            no_of_A = np.count_nonzero(mol.symbols == element)
            mol_energy -= (no_of_A*atom_energy[element])
        energy.append(mol_energy/Hartree)
    return np.array(energy)

def get_energies_perAtom(mols, atom_energy, energy_keyword = "energy"):
    """Function to get atomization energy per atom. This function is used for comparing results
    from kQEq, so energy is returned in eV. Dictionary of energies of atoms in vacuum is needed.    

    Parameters
    ----------
        mols : list
            list of ASE atoms with energies
        atom_energy : dictionary
            dictionary of energies of atoms in vacuum
        energy_keyword : string
            word used for energy extraction (energies has to be in eV) 

    Returns
    -------
        np.array
            Numpy array containing the atomization energies per atom in eV
    """
    energy = []
    for mol in mols:
        mol_energy = mol.info[energy_keyword]
        for element in atom_energy:
            no_of_A = np.count_nonzero(mol.symbols == element)
            mol_energy -= (no_of_A*atom_energy[element])
        mol_energy = mol_energy#/Hartree
        energy.append(mol_energy/len(mol))
    return energy

def get_whole_energies(mols, energy_keyword = "energy"):
    """Function to get energies of systems. Energy is returned in eV.    

    Parameters
    ----------
        mols : list
            list of ASE atoms with energies
        energy_keyword : string
            word used for energy extraction (energies has to be in eV) 

    Returns
    -------
        np.array
            Numpy array containing the energies in eV
    """
    energy = []
    for mol in mols:
        mol_energy = mol.info[energy_keyword]
        energy.append(mol_energy)
    return np.array(energy)#/Hartree

def get_whole_energies_perAtom(mols, energy_keyword = "energy"):
    """Function to get energies of systems per atom. Energy is returned in eV.    

    Parameters
    ----------
        mols : list
            list of ASE atoms with energies
        energy_keyword : string
            word used for energy extraction (energies has to be in eV) 

    Returns
    -------
        np.array
            Numpy array containing the energies per atom in eV
    """
    energy = []
    for mol in mols:
        mol_energy = mol.info[energy_keyword]
        energy.append(mol_energy/len(mol))
    return energy

def _get_R(tmpmol):
    """Helper function to get the R matrix of COM shifted coordinates

           R translates a vector of charges into the dipole vector (in atomic units)
           R*q -> dipole

    Parameters
    ----------
        tmpmol : obj
            ASE atoms objects with initial charges

    Returns
    -------
        np.array
            Numpy array containing the transformation matrix
    """
    mol = tmpmol.copy()
    Nat = len(mol)
    com = mol.get_center_of_mass()
    mol.translate(-com)
    R = mol.get_positions()
    #q = mol.get_initial_charges()
    return np.transpose(R/Bohr)

def _block_diag_rect(all_rs,dim1,dim2):
    """Helper function to construct blocked matrix from list of rectangular matrices.

    Parameters
    ----------
        all_rs : list
            List of rectangular matrices
        dim1: int
            First dimension of blocked matrix
        dim2: int
            Second dimension of blocked matrix

    Returns
    -------
        np.array
            Numpy array containing the blocked matrix
    """
    R_bar = np.zeros((dim2,dim1))
    counti = 0
    countj = 0

    for i,r in enumerate(all_rs):
        for ir in range(r.shape[0]):
            for jr in range(r.shape[1]):
                R_bar[ir+counti,jr+countj] = r[ir,jr]
        counti += r.shape[0]
        countj += r.shape[1]+1
    return R_bar

def _block_diag(all_as,dim):
    """Helper function to construct blocked matrix from list of square matrices.

    Parameters
    ----------
        all_as : list
            List of square matrices
        dim: int
            Dimension of blocked matrix

    Returns
    -------
        np.array
            Numpy array containing the blocked matrix
    """
    A_bar = np.zeros((dim,dim))
    counti = 0
    countj = 0

    for i,a in enumerate(all_as):
        for ia in range(a.shape[0]):
            for ja in range(a.shape[1]):
                A_bar[ia+counti,ja+countj] = a[ia,ja]
        counti += a.shape[0]
        countj += a.shape[1]
    return A_bar

def calc_gamma(a1,a2):
    return 1.0/np.sqrt(a1**2+a2**2)

################################################################################################################################
################################################## PLOTS #######################################################################
################################################################################################################################
import matplotlib.pyplot as plt

def plot_basics(kqeq, ref, title = None,myXlabel = None, myYlabel = None , preset="dipole", size_fig = (8,5), fig_dpi = 100,xlim = [-8,8], ylim = [-8,8], save = None):
    if preset == "dipole":
        xlabel = 'Reference Dipole Vector Elements'
        ylabel = 'kQEq Dipole Vector Elements'
    elif preset == "charges":
        xlabel = 'Reference charges'
        ylabel = 'kQEq charges'
    elif preset == "energy":
        xlabel = "Reference DFT Energy per atom / eV"
        ylabel = 'kQEq Energy per atom / eV'
    else:
        xlabel = myXlabel
        ylabel = myYlabel
    plt.figure(figsize=size_fig, dpi=fig_dpi)
    plt.scatter(ref,kqeq,alpha = 0.5)
    plt.plot(ref,ref,color="black")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.title(title,fontsize=26)
    plt.xlabel(xlabel,fontsize=22)
    plt.ylabel(ylabel,fontsize=22)
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    plt.show()


def plot_adv(kqeq, ref, title = None,myXlabel = None, myYlabel = None , preset="dipole", size_fig = (8,5), fig_dpi = 100,xlim = [-8,8], ylim = [-8,8], save = False):
    if preset == "dipole":
        xlabel = 'Reference Dipole Vector Elements'
        ylabel = 'kQEq Dipole Vector Elements'
    elif preset == "charges":
        xlabel = 'Reference charges'
        ylabel = 'kQEq charges'
    elif preset == "energy":
        xlabel = "Reference DFT Energy per atom / eV"
        ylabel = 'kQEq Energy per atom / eV'
    else:
        xlabel = myXlabel
        ylabel = myYlabel
    plt.figure(figsize=size_fig, dpi=fig_dpi)
    plt.scatter(ref,kqeq,alpha = 0.5)
    plt.plot([-10,10],[-10,10],color = "black")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.title(title,fontsize=26)
    plt.xlabel(xlabel,fontsize=22)
    plt.ylabel(ylabel,fontsize=22)
    #plt.savefig("trainDip.pdf",bbox_inches='tight')
    plt.show()












def create_lammps(mol, charges, enegs, elements, masses, name,hardness):
    input_file = open(f"{name}.in","w")
    input_file.write("#    for atom in mol:\n")
    input_file.write("units           metal\n")
    input_file.write("dimension       3\n")
    input_file.write("boundary        s s s\n")
    input_file.write("atom_style      charge\n")
    input_file.write("log log.lammps\n")
    input_file.write("fix enegID all property/atom d_eneg d_hard d_size\n")
    input_file.write(f"read_data {name}.data fix enegID NULL Enegs\n")
    input_file.write("pair_style      coul/gauss/cut 20.0\n")
    input_file.write("pair_coeff * *\n")
    #input_file.write(f"pair_coeff 1 1 20 {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['Zn'])}\n") # this for working long range
    #input_file.write(f"pair_coeff 2 2 20 {calc_gamma(uff_radius_qeq['O'],uff_radius_qeq['O'])}\n") # this for working long range
    #input_file.write(f"pair_coeff 1 2 20 {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['O'])}\n") # this for working long range
    input_file.write(f"pair_coeff 1 1 20\n")# {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['Zn'])}\n")
    input_file.write(f"pair_coeff 2 2 20\n")# {calc_gamma(uff_radius_qeq['O'],uff_radius_qeq['O'])}\n")
    input_file.write(f"pair_coeff 1 2 20\n")# {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['O'])}\n")
    #input_file.write(f"pair_coeff 1 1 20 {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['Zn'])}\n")
    #input_file.write(f"pair_coeff 2 2 20 {calc_gamma(uff_radius_qeq['O'],uff_radius_qeq['O'])}\n")
    #input_file.write(f"pair_coeff 1 2 20 {calc_gamma(uff_radius_qeq['Zn'],uff_radius_qeq['O'])}\n")
    input_file.write("dump 4a all custom 1 dump.myforce id type fx fy fz\n")
    input_file.write("dump myDump all atom 10 dump.lammpstrj\n")
    input_file.write("dump enegDump all custom 1 dump.eneg d_eneg\n")
    input_file.write("run 0\n")
    input_file.close()

    nAt = len(mol)
    lammps_elements = {}
    count_elements = 1
    for el in elements:
        lammps_elements[el] = count_elements
        count_elements += 1
    data_file = open(f"{name}.data","w")
    data_file.write("Comment lind\n")
    data_file.write("\n")
    data_file.write(f"          {nAt} atoms\n")
    data_file.write(f"          {len(elements)} atom types\n")
    data_file.write("\n")
    data_file.write("  -10.0   10.0      xlo xhi\n")
    data_file.write("  -10.0   10.0      ylo yhi\n")
    data_file.write("  -10.0   10.0      zlo zhi\n")
    data_file.write("\n")
    data_file.write(" Masses\n")
    data_file.write("\n")
    data_file.write("    1 65.38 \n    2 15.999\n")
    data_file.write("\n")
    data_file.write(" Atoms\n")
    data_file.write("\n")
    countAt = 1
    for at in mol:
        xyz = at.position
        symbol = lammps_elements[at.symbol]
        data_file.write(f"{countAt} {symbol} {charges[countAt-1]} {xyz[0]} {xyz[1]} {xyz[2]}\n")
        countAt += 1    
    data_file.write("\n")
    data_file.write(" Enegs\n")
    data_file.write("\n")
    countAt = 1
    for at in mol:
        data_file.write(f"{countAt} {enegs[countAt-1]} {hardness[at.symbol]} {uff_radius_qeq[at.symbol]/Bohr}\n")
        countAt += 1    
    data_file.close()

