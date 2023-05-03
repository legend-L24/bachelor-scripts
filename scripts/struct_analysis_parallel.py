import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import get_distances
from lmp_output_handler import lammps_dump_file_2_ase
import logging
import matplotlib.pyplot as plt

"""ToDo
- function to plot conductivity for different building block percentages
"""


def get_partial_rdf(atoms, rmax, nbins, center_idx=None, rest_idx=None):
    '''
    Parameters:
    -----------
    atoms:      ASE atoms object
    rmax:       float, maximum distance in \AA to perform rdf
    nbins:      int, number of bins
    center_idx: list, list of indices which to reference over for rdf
                (e.g. mask of all Al atoms)
    rest_idx:   list, list of indices which to account for in rdf
                (e.g. mask of all Li atoms)

    Returns:
    --------
    rdf:        np.array, with rdf values
    dists:      np.array, with distance values
    '''
    import numpy as np
    import math

    dm = atoms.get_all_distances(mic=True)
    rdf = np.zeros(nbins + 1)
    dr = float(rmax / nbins)

    center_idx = np.array(center_idx)
    r_ls = []
    for i in center_idx:
        for j in np.array(rest_idx):
            if i==j: continue
            
            rij = dm[i][j]
            r_ls.append(rij)
            #print("this is rij", i, j, rij)
            if rij != 0:
                index = int(math.ceil(rij / dr))
                if index <= nbins:
                    rdf[index] += 1
    #print("this is minimum of rij", min(r_ls))
    num=0
    for rij in r_ls:
        if 2>=rij>=0:
            num+=1
            print("wrong", rij)

    print("this is wrong possibility", num/len(r_ls))
    #print(r_ls)
    dists = []
    for i in range(1, nbins + 1):
        rrr = (i - 0.5) * dr
        dists.append(rrr)
        # Normalize with spherical approximation
        outer = 4 / 3 * np.pi * (rrr + dr) ** 3
        inner = 4 / 3 * np.pi * (rrr) ** 3
        vol = outer - inner
        rdf[i] /= vol

    return np.array(rdf[1:]), np.array(dists)


def get_element_atoms(atoms, el="P", idx=False):
    """get atoms object containing only atoms of type el"""
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
    """get 3x3x3 of atoms, handy to avoid problems with pbc"""
    new_atoms = Atoms(cell=atoms.cell, pbc=[1,1,1])
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                mult_factors = np.array([i,j,k])
                atoms_mult = atoms.copy()
                atoms_mult.positions+=np.dot(atoms.cell, mult_factors)
                new_atoms.extend(atoms_mult)
    return new_atoms


def get_building_blocks_1(atoms):
    """get building block percentages from atoms object"""
    comp = {"ps4":0, "p2s6":0, "p2s7":0}
    patoms = get_element_atoms(atoms)
    patoms2 = get_3xx3_cube(patoms)
    #print("check the code", len(patoms), len(patoms2))
    for p1 in patoms:
        for p2 in patoms2:
            dist = np.linalg.norm(p1.position-p2.position)
            #if dist < 2.5: print(dist)
            if 1. < dist < 2.5:
                comp["p2s6"]+=1
            elif 2.5 < dist < 4.5:
                comp["p2s7"]+=1
    comp["ps4"] = len(patoms) - comp["p2s6"] - comp["p2s7"]
    return [comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)]

def get_building_blocks(atoms):
    """get building block percentages from atoms object"""
    comp = {"ps4":0, "p2s6":0, "p2s7":0}
    patoms = get_element_atoms(atoms)
    patoms2 = get_3xx3_cube(patoms)
    #print("check the code", len(patoms), len(patoms2))
    for p1 in patoms:
        comp_p2s6, comp_p2s7 = 0,0
        for p2 in patoms2:
            dist = np.linalg.norm(p1.position-p2.position)
            #if 1<dist<2.5:print(dist)
            #if dist < 2.5: print(dist)
            #print("something happens")
            if 1. < dist < 2.3:    
                #print(dist)
                #print(dist)
                comp_p2s6 =1
            elif 2.3 < dist < 4.5:
                #print(dist)
                comp_p2s7 =1
        if comp_p2s7 and comp_p2s6:
            #print("one time meet comp_p2s6")
            comp["p2s6"]+=1; 
            #print("it is common");continue
        else:
            comp["p2s6"]+=comp_p2s6
            comp["p2s7"]+=comp_p2s7
    comp["ps4"] = len(patoms) - comp["p2s6"] - comp["p2s7"]
    return [comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)]

def get_building_blocks_S(atoms):
    """get building block percentages from atoms object"""
    comp = {"ps4":0, "p2s6":0, "p2s7":0}
    patoms = get_element_atoms(atoms, el='S')
    patoms2 = get_3xx3_cube(patoms)
    #print("check the code", len(patoms), len(patoms2))
    for p1 in patoms:
        comp_p2s6, comp_p2s7 = 0,0
        for p2 in patoms2:
            dist = np.linalg.norm(p1.position-p2.position)
            #if dist < 2.5: print(dist)
            if 1. < dist < 2.5:
                comp_p2s6 =1
            elif 2.5 < dist < 4.5:
                #print(dist)
                comp_p2s7 =1
        if comp_p2s7 and comp_p2s6:
            comp["p2s6"]+=1; 
            #print("it is common");continue
        else:
            comp["p2s6"]+=comp_p2s6
            comp["p2s7"]+=comp_p2s7
    comp["ps4"] = len(patoms) - comp["p2s6"] - comp["p2s7"]
    print([comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)])
    return [comp["ps4"]/len(patoms), comp["p2s6"]/len(patoms), comp["p2s7"]/len(patoms)]

def get_building_blocks_new(atoms):
    """get building block percentages from atoms object"""
    comp = {"ps3":0, "ps4":0, "p2s6_":0, "p2s7 or p2s6":0} #p2s6_ means the weired ion p2s6 4-
    patoms = get_element_atoms(atoms)
    satoms = get_element_atoms(atoms, el='S')
    pp_dist = patoms.get_all_distances()
    pp_dist += np.eye(pp_dist.shape[0], pp_dist.shape[1])*10
    ps_dist = get_distances(patoms.get_positions(), satoms.get_positions(), cell=patoms.cell, pbc=patoms.pbc)[-1]
    #print(ps_dist)
    #print("check the position, ps_dist shape, n(p), n(s): ",ps_dist.shape, len(patoms), len(satoms))
    ps_count = np.sum(ps_dist<2.6, axis=1)
    for i in range(len(patoms)):
        pp_min = np.min(pp_dist[i])
        if 1.0 < pp_min <= 2.5:
            comp["p2s6_"]+=1
            continue
        if 2.5 <= pp_min < 4.5:
            comp["p2s7 or p2s6"]+=1
            continue
        if 4.5 <= pp_min:
            if ps_count[i] == 3:comp["ps3"] += 1
            elif ps_count[i] == 4:comp["ps4"] += 1
            else: print("check weired number of S atoms :", ps_count[i])
            continue
        print("chekc one pair of P-P is too close, the distance: ", pp_min)
    #print("check the code", len(patoms), len(patoms2))
    return [comp["ps3"]/len(patoms), comp["ps4"]/len(patoms), comp["p2s6_"]/len(patoms), comp["p2s7 or p2s6"]/len(patoms)]


def read_dump(dumpfile):
    d_species = {2:15,1:11,3:16}; l_remove_type=[4]
    atoms = lammps_dump_file_2_ase(dumpfile,d_species,l_remove_type)
    return atoms
    #write(dumpfile[:-4]+"xyz", atoms)
'''
atoms = read("/home/yli/crystal_na3ps4/alpha_na3ps4.cif")
blocks = get_building_blocks_new(atoms)
print(blocks)
atoms = read("/home/yli/crystal_na3ps4/beta_na3ps4.cif")
blocks = get_building_blocks_new(atoms)
print(blocks)
atoms = read("/home/yli/crystal_na3ps4/gama_na3ps4.cif")
blocks = get_building_blocks_new(atoms)
print(blocks)
atoms = read("/home/yli/Full_Geo.xyz","::")
blocks = get_building_blocks_new(atoms[0])
print(blocks)
blocks = get_building_blocks(atoms[0])
print(blocks)
'''
from random import randint
import numpy as np
eval=False
if __name__=="__main__":
    '''
    blocks_ls = []
    atoms = read(f"/home/yli/na4p2s7_traj.xyz","::")
    print("this is atoms length", len(atoms))
    for atom in atoms:
        blocks = get_building_blocks_S(atom)
        blocks_ls.append([600, 0] + blocks)

        print(blocks)

    atoms = read(f"/home/yli/na3ps4_traj.xyz","::")
    print("this is atoms length", len(atoms))
    for atom in atoms:
        blocks = get_building_blocks_S(atom)
        blocks_ls.append([600, 0] + blocks)

        print(blocks)
    blocks_arr = np.matrix(blocks_ls)
    np.savetxt(f"/home/yli/database/glassy_na3ps4_rdf/67_75_comp.txt", blocks_arr)
    '''
    
    blocks_ls = []

    #for i in [0,1,3,5]:
    for i in [2,3,4,5]:
        #atoms = read(f"/home/yli/database/0320_comp//last/na4p2s7_1/{i}_comp.xyz","::")
        atoms = read(f"/home/yli/database/equ/na4p2s7/{i}_comp.xyz","300::")
        #atoms = read(f"/home/yli/database/0324/na4p2s7_1/na4p2s7_0318.xyz","::")
        print("this is atoms length", len(atoms))
        '''
        for atom in atoms:
            blocks = get_building_blocks_new(atom)
            frac = (randint(0,3)+10)/30
            blocks[-2]+=blocks[-1]*frac
            blocks[-1]-=blocks[-1]*frac
            blocks_ls.append([i*100+300, 0] + blocks[1:])
            #print(f"{i*50+350}K")
            print(blocks)
        '''
        for atom in atoms:
            blocks = get_building_blocks(atom)
            blocks_ls.append([i*100+300, 0] + blocks[:])
            #print(f"{i*50+350}K")
            print("this is one time:",blocks)
    blocks_arr = np.matrix(blocks_ls)
    print(blocks_arr.shape)
    np.savetxt(f"/home/yli/database/glassy_na3ps4_rdf/na4p2s7_0412_large.txt", blocks_arr)
        #print(np.max(blocks_arr, axis=0), np.min(blocks_arr, axis=0))
            #print(blocks)
    
    '''
    unit1_ls = [[],[],[],[],[],[]]
    unit2_ls = [[],[],[],[],[],[]]
    unit3_ls = [[],[],[],[],[],[]]
    for idx, ls in enumerate([unit1_ls, unit2_ls, unit3_ls]):
        for i in range(1,6):
            atoms = read(f"/home/yli/database/glassy_na3ps4_rdf/{i}.xyz","::")
            for atom in atoms:
                blocks = get_building_blocks(atom)
                ls[i].append(blocks[idx-1])
                #print(blocks)
    print(unit1_ls, unit2_ls, unit3_ls)
    '''
    '''
    plt.xlabel("Temperature / K")
    plt.ylabel("at.% P of species")
    plt.violinplot(unit1_ls[0])
    '''
            #print(blocks)
    #atoms = read_dump("/media/yli/LaCie/na3ps4/1/geom.dump")
    '''
    atoms = read("/home/yli/crystal_na3ps4/alpha_na3ps4.cif")
    blocks = get_building_blocks(atoms)
    print(blocks)
    atoms = read("/home/yli/crystal_na3ps4/beta_na3ps4.cif")
    blocks = get_building_blocks(atoms)
    print(blocks)
    atoms = read("/home/yli/crystal_na3ps4/gama_na3ps4.cif")
    blocks = get_building_blocks(atoms)
    print(blocks)
    '''
    '''
    atoms = read("/home/yli/na4p2s7.xyz","::")
    blocks = get_building_blocks(atoms[0])
    print(blocks)
    blocks = get_building_blocks(atoms[-1])
    print(blocks)
    blocks = get_building_blocks(atoms[10])
    print(blocks)
    '''
    
    #atoms = read("/media/yli/LaCie/final_compressed.xyz","::")
    nbins = 150
    
    atoms = [read("/home/yli/database/alpha_na3ps4.cif"),read("/home/yli/database/beta_na3ps4.cif"),read("/home/yli/database/gama_na3ps4.cif")]
    
    '''
    rdf_arr = np.zeros((nbins,))
    dist_arr = np.zeros((nbins,))
    for i in range(0, len(atoms)):
        
        p_ls = []
        s_ls = []
        na_ls = []
        for idx, ele in enumerate(atoms[0].get_chemical_symbols()):
            if ele=='P': p_ls.append(idx)
            if ele=='S': s_ls.append(idx)
            if ele=='Na': na_ls.append(idx)
        rdf_arr_0, dist_arr_0 = get_partial_rdf(atoms[0], rmax=6, nbins=nbins, center_idx=s_ls, rest_idx=p_ls)
        rdf_arr += rdf_arr_0; dist_arr = dist_arr_0

    plt.plot(dist_arr, rdf_arr/len(atoms), label="P-S")
    
    rdf_arr = np.zeros((nbins,))
    dist_arr = np.zeros((nbins,))
    for i in range(0, len(atoms)):
        
        p_ls = []
        s_ls = []
        na_ls = []
        for idx, ele in enumerate(atoms[0].get_chemical_symbols()):
            if ele=='P': p_ls.append(idx)
            if ele=='S': s_ls.append(idx)
            if ele=='Na': na_ls.append(idx)
        rdf_arr_0, dist_arr_0 = get_partial_rdf(atoms[0], rmax=6, nbins=nbins, center_idx=p_ls, rest_idx=p_ls)
        rdf_arr += rdf_arr_0; dist_arr = dist_arr_0

    plt.plot(dist_arr, rdf_arr/len(atoms), label="P-P")
    print(s_ls)
    '''
    nbins = 150
    '''
    rdf_arr = np.zeros((nbins,))
    dist_arr = np.zeros((nbins,))
    for i in range(0, len(atoms)):#len(atoms)):
        
        p_ls = []
        s_ls = []
        na_ls = []
        for idx, ele in enumerate(atoms[i].get_chemical_symbols()):
            if ele=='P': p_ls.append(idx)
            if ele=='S': s_ls.append(idx)
            if ele=='Na': na_ls.append(idx)
        rdf_arr_0, dist_arr_0 = get_partial_rdf(atoms[i], rmax=6, nbins=nbins, center_idx=get_element_atoms(atoms[i],'S',idx=True), rest_idx=get_element_atoms(atoms[i],'S',idx=True))
        rdf_arr += rdf_arr_0; dist_arr = dist_arr_0

    plt.plot(dist_arr, rdf_arr/len(atoms), label="S-S")

    plt.legend()
    plt.show()
    '''
    '''
    for i in range(0, len(atoms)-1):#len(atoms)):    
        rdf_arr_0, dist_arr_0 = get_partial_rdf(atoms[i], rmax=6, nbins=nbins, center_idx=get_element_atoms(atoms[i],'S',idx=True), rest_idx=get_element_atoms(atoms[i],'S',idx=True))
        plt.ylabel("RDF S-S"); plt.xlabel("r(A)")
        plt.plot(dist_arr_0, rdf_arr_0/len(atoms), label="S-S")
        plt.show()
        plt.clf()
    if eval:
        outfiles=["li3ps4_comp.txt","li7p3s11_comp.txt", "li4p2s7_comp.txt"]
        path = r"E:\MA_backup\tabea_ma_data\md_runs\production_runs\statistics_MD"
        crystals = ["Li3PS4", "Li7P3S11", "Li4P2S7"]
        Ts = [400,500, 600, 700,1000]
        folders = list(range(20))
        for crystal, outfile in zip(crystals, outfiles):
            with open(outfile, "w") as f:
                f.write("T ensemblecount ps4 p2s6 p2s7 \n")
            for T in Ts:
                for folder in folders:
                    atoms = read(r"{}\{}\{}\{}\geom.xyz".format(path, crystal, T, folder), ":")[-1]
                    blocks = get_building_blocks(atoms)
                    with open(outfile, "a") as f:
                        f.write("{} {} {} {} {}\n".format(T, folder, blocks[0], blocks[1], blocks[2]))
                        logging.info("wrote to file")
    '''
