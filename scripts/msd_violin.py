import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl
import sys

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
    lablist = [u"MSD(P)",u"P$_2$S$_6$$^{4-}$",u"P$_2$S$_7$$^{4-}$"]
    #pos=[400,500,600,700,800]
    pos=[500,600,700,800]#,450,550,650,750,850,950]
    #pos= ["5","6","7"]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,dataset in enumerate(datasets):
        print("this is dataset", dataset)
        print(f"the {i} th successes")
        ax.violinplot(dataset,pos,widths=50)
        #add_label(ax.violinplot(dataset,pos,widths=50), lablist[i])
    
    #ax.set_xticks([400,500,600],["one","two", "three"])
    #ax.set_xlim(400,700)
    #ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.set_xlabel("Temperature / K", fontsize = fs-3)
    ax.set_ylabel(r"MSD / $\rm \AA$^2", fontsize = fs-3)
    #ax.legend(*zip(*labels), bbox_to_anchor=(1.05,1.0), loc="upper left", prop = {"size": fs-5})
    #fig.subplots_adjust(right=0.8)
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


#path = sys.argv[1]+"/"
cmap = cm.get_cmap('coolwarm')
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({"figure.autolayout":True})
ratio_ls = []

p_msd = [[],[],[],[]]+[[],[],[],[],[],[]]
for phase in ["na3ps4","na3ps4_1","na3ps4_2"]:
#for phase in ["na3ps4_bridson","na3ps4_bridson_1","na3ps4_bridson_2"]:
    #phase = "na3ps4"#"na3ps4_bridson"
    '''
    files_p = [f"plots/{phase}_0_500K_p.txt",f"plots/{phase}_1_600K_p.txt",f"plots/{phase}_2_700K_p.txt",f"plots/{phase}_3_800K_p.txt",f"plots/{phase}_4_900K_p.txt"]
    files_s = [f"plots/{phase}_0_500K_s.txt",f"plots/{phase}_1_600K_s.txt",f"plots/{phase}_2_700K_s.txt",f"plots/{phase}_3_800K_s.txt",f"plots/{phase}_4_900K_s.txt"]
    files_li = [f"plots/{phase}_0_500K_li.txt",f"plots/{phase}_1_600K_li.txt",f"plots/{phase}_2_700K_li.txt",f"plots/{phase}_3_800K_li.txt",f"plots/{phase}_4_900K_li.txt"]
    '''
    files_p = [f"plots/{phase}_0_500K_p.txt",f"plots/{phase}_1_600K_p.txt",f"plots/{phase}_2_700K_p.txt",f"plots/{phase}_3_800K_p.txt",f"plots/{phase}_4_900K_p.txt"]
    files_s = [f"plots/{phase}_0_500K_s.txt",f"plots/{phase}_1_600K_s.txt",f"plots/{phase}_2_700K_s.txt",f"plots/{phase}_3_800K_s.txt",f"plots/{phase}_4_900K_s.txt"]
    files_li = [f"plots/{phase}_0_500K_li.txt",f"plots/{phase}_1_600K_li.txt",f"plots/{phase}_2_700K_li.txt",f"plots/{phase}_3_800K_li.txt",f"plots/{phase}_4_900K_li.txt"]
    dir = "/home/yli/database/"

    files_p = [dir+i for i in files_p]
    files_s = [dir+i for i in files_s]
    files_li = [dir+i for i in files_li]

    files_li = files_li[0:4]#+files_li[2:3]
    files_p = files_p[0:4]#+files_p[2:3]
    files_s = files_s[0:4]#+files_s[2:3]
    T_colorbar = np.linspace(500,900, 256)
    Ts = [500,600,700, 800]
    ratio = []
    #fig, ax = plt.subplots(figsize=(5,5))
    for i,file in enumerate(files_p):
        data = np.loadtxt(files_p[i]).T
        data_s = np.loadtxt(files_s[i]).T
        data_li = np.loadtxt(files_li[i]).T
        print(data.shape)
        p_msd[i].extend(list(data_s[690:700,1]))
        ratio.append(data[-1,1]/data_s[-1,1])
    ratio_ls.append(ratio)

for index, phase in enumerate(["na3ps4","na3ps4_1","na3ps4_2"]):
    #phase = "na3ps4"#"na3ps4_bridson"
    '''
    files_p = [f"plots/{phase}_0_500K_p.txt",f"plots/{phase}_1_600K_p.txt",f"plots/{phase}_2_700K_p.txt",f"plots/{phase}_3_800K_p.txt",f"plots/{phase}_4_900K_p.txt"]
    files_s = [f"plots/{phase}_0_500K_s.txt",f"plots/{phase}_1_600K_s.txt",f"plots/{phase}_2_700K_s.txt",f"plots/{phase}_3_800K_s.txt",f"plots/{phase}_4_900K_s.txt"]
    files_li = [f"plots/{phase}_0_500K_li.txt",f"plots/{phase}_1_600K_li.txt",f"plots/{phase}_2_700K_li.txt",f"plots/{phase}_3_800K_li.txt",f"plots/{phase}_4_900K_li.txt"]
    '''
    files_p = [f"plots/{phase}_0_450K_p.txt",f"plots/{phase}_1_550K_p.txt",f"plots/{phase}_2_650K_p.txt",f"plots/{phase}_3_750K_p.txt",f"plots/{phase}_4_850K_p.txt", f"plots/{phase}_5_950K_p.txt"]
    files_s = [f"plots/{phase}_0_450K_s.txt",f"plots/{phase}_1_550K_s.txt",f"plots/{phase}_2_650K_s.txt",f"plots/{phase}_3_750K_s.txt",f"plots/{phase}_4_850K_s.txt", f"plots/{phase}_5_950K_p.txt"]
    files_li = [f"plots/{phase}_0_450K_li.txt",f"plots/{phase}_1_550K_li.txt",f"plots/{phase}_2_650K_li.txt",f"plots/{phase}_3_750K_li.txt",f"plots/{phase}_4_850K_li.txt", f"plots/{phase}_5_950K_li.txt"]
    dir = "/home/yli/database/"

    files_p = [dir+i for i in files_p]
    files_s = [dir+i for i in files_s]
    files_li = [dir+i for i in files_li]

    files_li = files_li[0:6]#+files_li[2:3]
    files_p = files_p[0:6]#+files_p[2:3]
    files_s = files_s[0:6]#+files_s[2:3]
    T_colorbar = np.linspace(500,900, 256)
    Ts = [450,550,650,750,850,950]
    ratio = []
    #fig, ax = plt.subplots(figsize=(5,5))
    for i,file in enumerate(files_p):
        data = np.loadtxt(files_p[i]).T
        data_s = np.loadtxt(files_s[i]).T
        data_li = np.loadtxt(files_li[i]).T
        #idx = i+4
        #print(f"this is {i}th:",data[-5:-1,1])
        print(data.shape)
        p_msd[i+4].extend(list(data_li[600:610,1]))
        ratio.append(data[-1,1]/data_s[-1,1])
    ratio_ls[index]+=ratio

plt.title(phase)
plt.xlabel("T / K")
plt.ylabel("MSD(P)/MSD(S)")
ratio_ls = [[i[4],i[0],i[5],i[1],i[6],i[2],i[7],i[3],i[8],i[9]] for i in ratio_ls]
for ratio in ratio_ls:
    plt.plot([450,500,550,600,650,700,750,800,850,950],ratio)
plt.show()

labels = []
print(len(p_msd), len(p_msd[0]))
#p_msd_arr

do_violin_plot([p_msd], savename="{}{}".format(file[:-3],".pdf"))
plt.title("MSD of S")
plt.xticks(np.arange(400, 1000, 50), fontsize=15)
plt.show()
'''
plt.xlabel("T / K")
plt.ylabel("MSD(S)/MSD(Na)")
for ratio in ratio_ls:
    plt.plot(Ts,ratio)
'''
#plt.show()
#plt.savefig(path+phase+"_new_Na.png")