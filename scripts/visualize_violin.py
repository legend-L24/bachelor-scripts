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


def add_label(violin, label):
    print(violin["bodies"][0].get_facecolor())
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))
    return color, labels

def add_label_1(violin, label=None, color=None):
    if label:
        labels.append((mpatches.Patch(color=color), label))
    else:
        labels.append((mpatches.Patch(color=color)))
def do_violin_plot(datasets, savename="violinplot.pdf"):
    fs = 18
    lablist = [u"PS$_4$$^{3-}$",u"P$_2$S$_6$$^{4-}$",u"P$_2$S$_7$$^{4-}$"]
    #pos=[400,500,600,700,800]
    pos=[500,600,700,800]
    #pos= ["5","6","7"]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,dataset in enumerate(datasets):
        add_label(ax.violinplot(dataset,pos,widths=50), lablist[i])
    
    #ax.set_xticks([400,500,600],["one","two", "three"])
    #ax.set_xlim(400,700)
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.set_xlabel("Temperature / K", fontsize = fs-3)
    ax.set_ylabel(r"At % p of species", fontsize = fs-3)
    ax.legend(*zip(*labels), bbox_to_anchor=(1.05,1.0), loc="upper left", prop = {"size": fs-5})
    fig.subplots_adjust(right=0.8)
    plt.savefig(savename)
def do_violin_plot_1D(datasets, savename="violinplot.pdf"):
    fs = 18
    lablist = ["Initial_density","Final_density"]
    #pos=[400,500,600,700,800]
    pos=[400,500,600,700]
    #pos= ["5","6","7"]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,dataset in enumerate(datasets):
        add_label(ax.violinplot(dataset,pos,widths=50), lablist[i])
    
    #ax.set_xticks([400,500,600],["one","two", "three"])
    #ax.set_xlim(400,700)
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.set_xlabel("Temperature / K", fontsize = fs-3)
    ax.set_ylabel(r"Volume/ A^3", fontsize = fs-3)
    ax.legend(*zip(*labels), loc="upper left", prop = {"size": fs-4})
    plt.savefig(savename)

def do_violin_plot_compared(datasets, savename="violinplot.pdf"):
    fs = 18
    lablist = [u"PS3$^{-}$",u"PS4$^{3-}$",u"P2S6$^{4-}$",u"P2S7$^{4-}$ and P2S6 $^{2-}$"]
    #pos=[400,500,600,700,800]
    pos= [500, 600, 700, 800, 900, 1000]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color_ls = []
    #color_ls = ['g','b','y','b']
    for i,dataset in enumerate(datasets[0:4]):
        color_ls.append(add_label(ax.violinplot(dataset,pos[:2], widths=50), lablist[i]))
    print("check color", color_ls)
    ax.set_xlabel("Temperature / K", fontsize = fs)
    ax.set_ylabel(r"At % p of species", fontsize = fs)
    ax.legend(*zip(*labels), loc="right", prop = {"size": fs-4})
    for i,dataset in enumerate(datasets[4:]):
        print(dataset.shape, len(pos[2:]), i)
        add_label_1(ax.violinplot(dataset,pos[2:], widths=50), color=color_ls[i])
    x_label = [u"Na$_{4}$P$_{2}$S$_{7}$-GAP", u"Na$_{3}$PS$_{4}$-GAP", u"Na$_{4}$P$_{2}$S$_{7}$-AIMD", u"Na$_{3}$PS$_{4}$-AIMD",u"Na$_{4}$P$_{2}$S$_{7}$-NMR",  u"Na$_{3}$PS$_{4}$-NMR"]
    ax.set_xticks(pos,x_label, fontsize=fs-4)
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    #ax.legend()
    plt.tight_layout()
    plt.savefig(savename)
def do_violin_plot_compared_new(datasets, savename="violinplot.pdf"):
    fs = 18
    lablist = [u"PS3$^{-}$",u"PS4$^{3-}$",u"P2S6$^{4-}$",u"P2S7$^{4-}$ and P2S6 $^{2-}$"]
    #pos=[400,500,600,700,800]
    pos= [500, 600, 700, 800, 900, 1000]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color_ls = []
    #color_ls = ['g','b','y','b']
    for i,dataset in enumerate(datasets):
        print(dataset.shape)
        color_ls.append(add_label(ax.violinplot(dataset,pos, widths=50), lablist[i]))
    print("check color", color_ls)
    ax.set_xlabel("Temperature / K", fontsize = fs)
    ax.set_ylabel(r"At % p of species", fontsize = fs)
    ax.legend(*zip(*labels), loc="right", prop = {"size": fs-4})
    x_label = [u"Na$_{4}$P$_{2}$S$_{7}$-GAP", u"Na$_{3}$PS$_{4}$-GAP", u"Na$_{4}$P$_{2}$S$_{7}$-AIMD", u"Na$_{3}$PS$_{4}$-AIMD",u"Na$_{4}$P$_{2}$S$_{7}$-NMR",  u"Na$_{3}$PS$_{4}$-NMR"]
    ax.set_xticks(pos,x_label, fontsize=fs-4)
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    #ax.legend()
    plt.tight_layout()
    plt.savefig(savename)
'''
for file in ["li3ps4_comp.txt", "li7p3s11_comp.txt", "li4p2s7_comp.txt"]:
    labels = []
    data= np.loadtxt(file, skiprows=1)[:,2:]
    ps4 =  data[:,0].reshape((5, int(len(data)/5))).T
    p2s6 = data[:,1].reshape((5, int(len(data)/5))).T
    p2s7 = data[:,2].reshape((5, int(len(data)/5))).T
    do_violin_plot([ps4,p2s6,p2s7], savename="{}{}".format(file[:-3],".pdf"))
plt.show()
'''
'''
for file in ["/home/yli/database/glassy_na3ps4_rdf/67_75_comp.txt"]:
    labels = []
    data= np.loadtxt(file, skiprows=0)[:,2:]
    print(data.shape, len(data))
    #na4p2s7 = data[:1+len(data)/2,:].reshape(4, len(data)/2).T
    #na3ps4 = data[:1+len(data)/2,:].reshape(4, len(data)/2).T
    ps3 = data[:,0].reshape(2, int(len(data)/2)).T
    ps4 = data[:,1].reshape(2, int(len(data)/2)).T
    p2s6 = data[:,2].reshape(2, int(len(data)/2)).T
    p2s7 = data[:,3].reshape(2, int(len(data)/2)).T
    #aimd_na4p2s7 = np.array([0.05, 0.22, 0.09, 0.58]).reshape(4,1)
    #aimd_na3ps4 = np.array([5.5, 0.67, 0.11, 0.11]).reshape(4,1)
    #exp_na4p2s7 = np.array([0, 0.09, 0.19, 0.68]).reshape(4,1)
    #exp_na3ps4 = np.array([0, 0.70, 0.09, 0.18]).reshape(4,1)
    ps3_1 = np.array([0.05, 0.055, 0.0, 0]*20).reshape(20,4)
    ps4_1 = np.array([0.22, 0.67, 0.09, 0.70]*20).reshape(20,4)
    p2s6_1 = np.array([0.09, 0.105, 0.19, 0.09]*20).reshape(20,4)
    p2s7_1 = np.array([0.58, 0.115, 0.68, 0.18]*20).reshape(20,4)
    #do_violin_plot_compared_new([ps3, ps4, p2s6, p2s7, ps3_1, ps4_1, p2s6_1, p2s7_1], savename="{}{}".format(file[:-3],".pdf"))
    do_violin_plot_compared_new([np.hstack((ps3,ps3_1)), np.hstack((ps4,ps4_1)), np.hstack((p2s6,p2s6_1)), np.hstack((p2s7,p2s7_1))], savename="{}{}".format(file[:-3],".pdf"))

plt.show()
'''
'''
for file in ["/home/yli/database/glassy_na3ps4_rdf/na4p2s7_2_comp.txt"]:
    labels = []
    data= np.loadtxt(file, skiprows=0)[:,2:]
    print(data.shape)
    ps3 =  data[:,0].reshape((3, int(len(data)/3))).T
    ps4 =  data[:,1].reshape((3, int(len(data)/3))).T
    p2s6 = data[:,2].reshape((3, int(len(data)/3))).T
    p2s7 = data[:,3].reshape((3, int(len(data)/3))).T
    do_violin_plot([ps3, ps4,p2s6,p2s7], savename="{}{}".format(file[:-3],".pdf"))
plt.show()
'''

# The violin plot for volume change

if 0:
    labels = []
    vol_ls_initial = [5959.52301645441, 5970.111725440384, 5794.457160928438, 5959.52301645441, 5970.111725440384, 5794.457160928438, 5959.52301645441, 5970.111725440384, 5794.457160928438, 5959.52301645441, 5970.111725440384, 5794.457160928438]
    vol_ls_final = [5890.752312064266, 5925.92084965482, 5809.459410402131, 5862.415893905016, 5982.808621257674, 5889.833155668627, 5906.578117389765, 5998.522301842431, 5972.762806836655, 5903.18176956816, 6149.188352125968, 6118.03630570104]

    vol_ls_initial = [6204.084980746079, 6211.559802620195, 6207.482880409858, 6204.084980746079, 6211.559802620195, 6207.482880409858, 6204.084980746079, 6211.559802620195, 6207.482880409858, 6204.084980746079, 6211.559802620195, 6207.482880409858]
    vol_ls_final = [6273.158438386751, 6558.398445343162, 6372.883987702624, 6435.304789686414, 6423.770700154432, 6556.062998823564, 6890.288853268668, 6737.114739113809, 6781.190144102534, 6992.316012728252, 6960.274176794933, 7120.846070921048]
    data= np.array(vol_ls_initial)
    print(data.shape)
    vol_arr_0 =  data.reshape((4, int(len(data)/4))).T
    data= np.array(vol_ls_final)
    print(data.shape)
    vol_arr_1 =  data.reshape((4, int(len(data)/4))).T
    do_violin_plot_1D([vol_arr_0, vol_arr_1])
    plt.show()

# The curve show how the compositions change with time every 50 ps lablist = [u"PS$_4$$^{3-}$",u"P$_2$S$_6$$^{4-}$",u"P$_2$S$_7$$^{4-}$"]
if 0:
    fs = 32
    lw2 = 2
    plt.figure(figsize=[38, 8], dpi=100)
    for file in ["/home/yli/database/glassy_na3ps4_rdf/na3ps4_0324_simple.txt"]:
        labels = []
        data= np.loadtxt(file, skiprows=0)[:,2:]
        
        print(data.shape)
        ps4 =  data[:,0].reshape((4, int(len(data)/4))).T*100
        p2s6 = data[:,1].reshape((4, int(len(data)/4))).T*100
        p2s7 = data[:,2].reshape((4, int(len(data)/4))).T*100
        print(ps4.shape)
        T_ls = [500, 600, 700, 800]
        for idx, T in enumerate(T_ls):
            ax = plt.subplot(1,len(T_ls),idx+1)
            ax.spines['bottom'].set_linewidth(lw2)
            ax.spines['left'].set_linewidth(lw2)
            ax.spines['right'].set_linewidth(lw2)
            ax.spines['top'].set_linewidth(lw2)
            num = p2s7.shape[0]
            plt.title(u"Composition for Na$_4$P$_2$S$_{7}$"+f" in {T}K")
            plt.xlabel("Simulation time/ps", fontsize = fs-10)
            plt.ylabel(r"At % p of species", fontsize = fs-10)
            plt.plot(np.arange(0, num)*50, ps4[:,idx], label=u"PS$_4$$^{3-}$")
            plt.plot(np.arange(0, num)*50, p2s6[:,idx], label=u"P$_2$S$_6$$^{4-}$")
            plt.plot(np.arange(0, num)*50, p2s7[:,idx], label=u"P$_2$S$_7$$^{4-}$")
            plt.legend(prop={"size":fs-12},)
        plt.show()


for file in ["/home/yli/database/glassy_na3ps4_rdf/na3ps4_bridson_final.txt"]:
    labels = []
    data= np.loadtxt(file, skiprows=0)[:80,2:]
    print(data.shape)
    ps4 =  data[:,0].reshape((4, int(len(data)/4))).T*100
    p2s6 = data[:,1].reshape((4, int(len(data)/4))).T*100
    p2s7 = data[:,2].reshape((4, int(len(data)/4))).T*100
    print(ps4)
    #do_violin_plot([ps4,p2s6,p2s7], savename="{}{}".format(file[:-3],".pdf"))
    do_violin_plot([ps4,p2s6,p2s7], savename="{}{}".format(file[:-3],".pdf"))
plt.show()
#plt.savefig("/home/yli/database/na3ps4_0307.png")
#plt.show()
