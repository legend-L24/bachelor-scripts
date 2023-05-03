from ase import Atoms
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams.update({"font.size":15})


def get_element_atoms(atoms, el="P", idx=False):
    """get """
    elatoms = Atoms(cell=atoms.cell, pbc=[1,1,1])
    idx_list =[]
    for i,atom in enumerate(atoms):
        if atom.symbol==el:
            elatoms.append(atom)
            idx_list.append(i)
    if idx:
        return idx_list
    else:
        return elatoms


def get_3xx3_cube(atoms):
    new_atoms = Atoms(cell=atoms.cell, pbc=[1,1,1])
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                mult_factors = np.array([i,j,k])
                atoms_mult = atoms.copy()
                atoms_mult.positions+=np.dot(atoms.cell, mult_factors)
                new_atoms.extend(atoms_mult)
    return new_atoms


def get_building_blocks(atoms):
    comp = {"ps4":0, "p2s6":0, "p2s7":0}
    patoms = get_element_atoms(atoms)
    patoms2 = get_3xx3_cube(patoms)
    for p1 in patoms:
        for p2 in patoms2:
            dist = np.linalg.norm(p1.position-p2.position)
            if 1. < dist < 2.5:
                comp["p2s6"]+=1
            elif 2.5 < dist < 4.5:
                comp["p2s7"]+=1
    comp["ps4"] = len(patoms) - comp["p2s6"] - comp["p2s7"]
    return [comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)]

def get_building_blocks_new(atoms):
    comp = {"ps3":0, "ps4":0, "p2s6":0, "p2s6_":0,"p2s7":0} # p2s6_ stand for the weired ions: p2s6 4-
    patoms = get_element_atoms(atoms)
    patoms2 = get_3xx3_cube(patoms)
    for p1 in patoms:
        for p2 in patoms2:
            dist = np.linalg.norm(p1.position-p2.position)
            if 1. < dist < 2.5:
                comp["p2s6"]+=1
            elif 2.5 < dist < 4.5:
                comp["p2s7"]+=1
    comp["ps4"] = len(patoms) - comp["p2s6"] - comp["p2s7"]
    return [comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)]



def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))


def do_violin_plot(datasets, savename="violinplot.pdf"):
    lablist = ["PS4","P2S6","P2S7"]
    pos=[400,500,600,700,1000]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,dataset in enumerate(datasets):
        add_label(ax.violinplot(dataset,pos, widths=50), lablist[i])
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.set_xlabel("Temperature / K")
    ax.set_ylabel(r"At % p of species")
    ax.legend(*zip(*labels), loc=2)
    plt.savefig(savename)


for file in ["li3ps4_comp.txt", "li7p3s11_comp.txt", "li4p2s7_comp.txt"]:
    labels = []
    data= np.loadtxt(file, skiprows=1)[:,2:]
    ps4 =  data[:,0].reshape((5, int(len(data)/5))).T
    p2s6 = data[:,1].reshape((5, int(len(data)/5))).T
    p2s7 = data[:,2].reshape((5, int(len(data)/5))).T
    do_violin_plot([ps4,p2s6,p2s7], savename="{}{}".format(file[:-3],".pdf"))
plt.show()
