import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams["legend.handlelength"] = .5
"""500K"""
#iface1 = read("plots/iface500iface1._")
#iface500crystal1._2__s
#iface2 = read()
#crystal = read()


def listplot(list):
    for el in list:
        plt.plot(el[:,0], el[:,1])
    plt.show()


def get_shortest_el(list):
    l = 100000
    for i,el in enumerate(list):
        if len(el) < l:
            l = len(el)
    return l


"""800K"""
space_resolved = False
if space_resolved:
    iface1 = np.loadtxt("../plots/ifaceiface1_0___li.txt").T
    iface2 = np.loadtxt("../plots/ifaceiface2_0___li.txt").T
    iface3 = np.loadtxt("../plots/ifaceiface1_1___li.txt").T
    iface4 = np.loadtxt("../plots/ifaceiface2_1___li.txt").T
    iface5 = np.loadtxt("../plots/ifaceiface1_2___li.txt").T
    iface6 = np.loadtxt("../plots/ifaceiface2_2___li.txt").T
    ifacelist = [iface1, iface2, iface3, iface4, iface5, iface6]
    listplot(ifacelist)
    l = get_shortest_el(ifacelist)
    iface_mat = np.array([iface1[:l,1],iface2[:l,1], iface3[:l,1], iface4[:l,1], iface5[:l,1], iface6[:l,1]])
    iface_av = np.mean(iface_mat, axis=0)
    iface_std = np.std(iface_mat, axis=0)

    cr1 = np.loadtxt("../plots/ifacecrystal1_0___li.txt").T
    cr2 = np.loadtxt("../plots/ifacecrystal1_1___li.txt").T
    cr3 = np.loadtxt("../plots/ifacecrystal1_2___li.txt").T
    crlist = [cr1, cr2, cr3]
    listplot(crlist)
    l = get_shortest_el(crlist)
    cr_mat = np.array([cr1[:l,1],cr2[:l,1], cr3[:l,1]])
    cr_av = np.mean(cr_mat, axis=0)
    cr_std = np.std(cr_mat, axis=0)

    gl1 = np.loadtxt("../plots/ifaceglass_0___li.txt").T
    gl2 = np.loadtxt("../plots/ifaceglass_1___li.txt").T
    gl3 = np.loadtxt("../plots/ifaceglass_2___li.txt").T
    gllist = [gl1, gl2, gl3]
    listplot(gllist)
    l = get_shortest_el(gllist)
    gl_mat = np.array([gl1[:l,1],gl2[:l,1], gl3[:l,1]])
    gl_av = np.mean(gl_mat, axis=0)
    gl_std = np.std(gl_mat, axis=0)

    #plt.errorbar(iface1[:l,0], iface_av, iface_std, label="interface")
    plt.errorbar(cr1[:l,0], cr_av, cr_std, label="crystal")
    plt.errorbar(gl1[:l,0], gl_av, gl_std, label="glass")
    plt.legend()
    plt.show()


#"""whole iface"""
#fig = plt.figure(figsize=(5,5))
#ax = fig.add_subplot(111)
#colors= ["blue", "green", "red"]
#for i in range(3):
#    for j in range(3):
#        data = np.loadtxt("../plots/iface/iface_compl_500__"+str(i)+"_"+str(j)+"__li.txt").T
#        if not (i==2 and j ==0):
#            plt.plot(data[:,0], data[:,1], color=colors[i])
#    plt.xlabel("t / ps")
#    plt.ylabel("msd / $\AA^2$")
#    ax.set_aspect(1.0 / ax.get_data_ratio())
#plt.show()


"""spatially resolved"""
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
labels = ["glass", "interface",  "crystal"]
colors =  ["mediumturquoise", "royalblue", "midnightblue"]

data = np.loadtxt("../plots/Li3PS4_16_500_glass__li.txt").T
plt.plot(data[:200,0], data[:200,1], color="mediumturquoise", label="amorphous Li$_3$PS$_4$", linewidth="5")
for i in range(3):
    for j in range(3):
        for k, type in enumerate(["glass.xyz", "iface1.xyz", "crystal1.xyz"]):
            if not (i == 2 and j == 0):
                data = np.loadtxt("../plots/iface/iface_"+type+"_"+str(i)+"_"+str(j)+"__li.txt").T
                if i==0 and j==0:
                    plt.plot(data[:, 0], data[:, 1], color=colors[k], label=labels[k])
                else:
                    plt.plot(data[:, 0], data[:, 1], color=colors[k])
data = np.loadtxt("../plots/beta_corr_lowT_500_li.txt").T
plt.plot(data[:200,0], data[:200,1], color="midnightblue", label="$\\beta$-Li$_3$PS$_4$", linewidth="5")
plt.legend()
plt.xlabel("t / ps")
plt.ylabel("msd / $\AA^2$")
ax.set_aspect(1.0 / ax.get_data_ratio())
plt.savefig("/home/huss/MA/Report/msd_iface.pdf")
plt.show()

