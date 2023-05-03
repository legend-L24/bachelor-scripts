#!/usr/bin/env python
# coding: utf-8

# ### Compare the two different sampling approaches (sampling on grid + voronoi and bridson sampling + fps) with respect to rdf and energy

# In[17]:


from partial_rdf import *


# In[18]:


"""load files"""
from ase.io import read, write
import quippy
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

voro_structs = read("../tests/final_compressed_voro.xyz", ":")
bridson_fps_structs = read("../tests/final_compressed_fps.xyz", ":")


# In[19]:


def get_density(atoms):
    """Get density of cell in g/cm^3"""
    masses = atoms.get_masses()
    N_A = 6.022 * 10 ** 23
    mass = sum([mass / N_A for mass in masses])
    volume = atoms.get_volume() * 10 ** (-24)
    density = mass / volume
    return(density)


# In[20]:


"""rdfs"""
for i, el1 in enumerate(["Na", "P", "S"]):
    for j, el2 in enumerate(["Na", "P", "S"]):
        if not j > i:
            rdf_av_fps, rdf_av_voro = [],[]
            for bridson_fps_struct, voro_struct in zip(bridson_fps_structs, voro_structs):
                rdf_fps, dist = get_partial_rdf(bridson_fps_struct, 6, 50, el1, el2)
                rdf_voro, dist = get_partial_rdf(voro_struct, 6, 50, el1, el2)
                rdf_av_fps.append(rdf_fps)
                rdf_av_voro.append(rdf_voro)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            rdf_av_fps = np.average(np.array(rdf_av_fps),axis=0)
            rdf_av_voro = np.average(np.array(rdf_av_voro), axis=0)
            plt.plot(dist, rdf_av_voro, label="grid voro")
            plt.plot(dist, rdf_av_fps, label="bridson fps")
            plt.xlabel("distance / $\AA$")
            plt.ylabel("rdf")
            plt.title(f"{el1}-{el2}")
            plt.legend()
            plt.tight_layout()
            plt.show()


# In[21]:


import quippy
"""energetics"""
pot = quippy.potential.Potential(param_filename="../tests/amorphous_lps_sampling_copy/gap/gp_2b_soap.xml")
fps_energies = []
voro_energies = []
for bridson_fps_struct in bridson_fps_structs:
    pot.calculate(bridson_fps_struct, ["forces", "energy"])
    fps_energies.append(pot.results["energy"]/len(bridson_fps_struct))
for voro_struct in voro_structs:
    pot.calculate(voro_struct, ["forces", "energy"])
    voro_energies.append(pot.results["energy"]/len(voro_struct))
plt.violinplot([fps_energies, voro_energies], positions=[0,1])
plt.xticks([0,1], ["bridson fps", "grid voro"], rotation=60)
plt.ylabel("energy per atom / eV")
plt.show()



# ### Compare after sampling and add computational sintering as possible "sampling" method

# In[1]:


from ase.io import read


def convert_symbols(str):
    dict = {"Li":"S", "He":"P", "H":"Na"}
    for key in dict.keys():
        str = str.replace(key, dict[key])
    return str

def read_dump(path, s=":"):
    atomslist = read(path, s)
    if not type(atomslist) is list:
        atomslist.symbols = convert_symbols(str(atomslist.symbols))
        return atomslist
    for i in range(len(atomslist)):
        atomslist[i].symbols = convert_symbols(str(atomslist[i].symbols))
    return atomslist



# In[23]:


"""energetics"""
import quippy
import matplotlib.pyplot as plt
pot = quippy.potential.Potential(param_filename="../tests/amorphous_lps_sampling_copy/gap/gp_2b_soap.xml")
fps_energies = []
voro_energies = []
melt_quench_energies = []
"""bridson"""
for i in range(3):
    struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/bridson_fps/{i}/geom.dump", "-1")
    pot.calculate(struct, ["forces", "energy"])
    fps_energies.append(pot.results["energy"]/len(struct))
"""grid/voro"""
for i in range(3):
    struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-1")
    pot.calculate(struct, ["forces", "energy"])
    voro_energies.append(pot.results["energy"]/len(struct))
"""melt-quench"""
struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-1")
pot.calculate(struct, ["forces", "energy"])
melt_quench_energies.append(pot.results["energy"]/len(struct))
# ToDo
plt.violinplot([fps_energies, voro_energies, melt_quench_energies], positions=[0,1,2])
plt.xticks([0, 1, 2], ["bridson fps", "grid voro", "melt-quench"], rotation=60)
plt.ylabel("energy per atom / eV")
plt.show()


# In[24]:


"""rdfs"""
# define voro structs
voro_structs = []
for i in range(3):
    voro_structs.extend(read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-5:"))
# define bridson structs
bridson_structs = []
for i in range(3):
    bridson_structs.extend(read_dump(f"/work/home/thuss/calc/naps/sampl_compar/bridson_fps/{i}/geom.dump", "-5:"))
#define melt-quench
melt_quench = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-1")

for i, el1 in enumerate(["Na", "P", "S"]):
    for j, el2 in enumerate(["Na", "P", "S"]):
        if not j > i:
            rdf_av_fps, rdf_av_voro = [],[]
            for bridson_struct, voro_struct in zip(bridson_structs, voro_structs):
                rdf_fps, dist = get_partial_rdf(bridson_struct, 6, 50, el1, el2)
                rdf_voro, dist = get_partial_rdf(voro_struct, 6, 50, el1, el2)
                rdf_av_fps.append(rdf_fps)
                rdf_av_voro.append(rdf_voro)
            rdf_melt, dist = get_partial_rdf(melt_quench, 6,50, el1, el2)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            rdf_av_fps = np.average(np.array(rdf_av_fps),axis=0)
            rdf_av_voro = np.average(np.array(rdf_av_voro), axis=0)
            plt.plot(dist, rdf_av_voro, label="grid voro")
            plt.plot(dist, rdf_av_fps, label="bridson fps")
            plt.plot(dist, rdf_melt, label="melt-quench")
            plt.xlabel("distance / $\AA$")
            plt.ylabel("rdf")
            plt.title(f"{el1}-{el2}")
            plt.legend()
            plt.tight_layout()
            plt.show()


# In[21]:


"""microchemistries"""
from microchemistry import get_building_blocks, get_el_atoms
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))

fps_bb = []
voro_bb = []
melt_quench_bb = []
"""bridson"""
for i in range(3):
    struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/bridson_fps/{i}/geom.dump", "-1")
    print(struct)
    fps_bb.append(get_building_blocks(struct))
fps_bb = np.array(fps_bb)
"""grid/voro"""
for i in range(3):
    struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-1")
    voro_bb.append(get_building_blocks(struct))
voro_bb = np.array(voro_bb)
"""melt-quench"""
struct = read_dump(f"/work/home/thuss/calc/naps/sampl_compar/grid_voro/{i}/geom.dump", "-1")
melt_quench_bb.append(get_building_blocks(struct))
melt_quench_bb = np.array(melt_quench_bb)
print(melt_quench_bb)
labels = []
add_label(plt.violinplot([fps_bb[:,0], voro_bb[:,0], melt_quench_bb[:,0]], positions=[0,1,2]), label="PS4")
add_label(plt.violinplot([fps_bb[:,1], voro_bb[:,1], melt_quench_bb[:,1]], positions=[0,1,2]), label="P2S6")
add_label(plt.violinplot([fps_bb[:,2], voro_bb[:,2], melt_quench_bb[:,2]], positions=[0,1,2]), label="P2S7")
plt.legend(*zip(*labels), loc=5)
plt.xticks([0, 1, 2], ["bridson fps", "grid voro", "melt-quench"], rotation=60)
plt.ylabel("building block occurence / %")
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




