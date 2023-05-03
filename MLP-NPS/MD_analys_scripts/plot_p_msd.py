import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl
import sys
path = sys.argv[1]+"/"
cmap = cm.get_cmap('coolwarm')
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({"figure.autolayout":True})

phase = "gama_new"#"na3ps4_bridson"
'''
files_p = [f"plots/{phase}_0_500K_p.txt",f"plots/{phase}_1_600K_p.txt",f"plots/{phase}_2_700K_p.txt",f"plots/{phase}_3_800K_p.txt",f"plots/{phase}_4_900K_p.txt"]
files_s = [f"plots/{phase}_0_500K_s.txt",f"plots/{phase}_1_600K_s.txt",f"plots/{phase}_2_700K_s.txt",f"plots/{phase}_3_800K_s.txt",f"plots/{phase}_4_900K_s.txt"]
files_li = [f"plots/{phase}_0_500K_li.txt",f"plots/{phase}_1_600K_li.txt",f"plots/{phase}_2_700K_li.txt",f"plots/{phase}_3_800K_li.txt",f"plots/{phase}_4_900K_li.txt"]
files_p = [f"plots/{phase}_0_400K_p.txt",f"plots/{phase}_1_600K_p.txt",f"plots/{phase}_2_800K_p.txt",f"plots/{phase}_3_1000K_p.txt",f"plots/{phase}_4_1200K_p.txt"]
files_s = [f"plots/{phase}_0_400K_s.txt",f"plots/{phase}_1_600K_s.txt",f"plots/{phase}_2_800K_s.txt",f"plots/{phase}_3_1000K_s.txt",f"plots/{phase}_4_1200K_s.txt"]
files_li = [f"plots/{phase}_0_400K_li.txt",f"plots/{phase}_1_600K_li.txt",f"plots/{phase}_2_800K_li.txt",f"plots/{phase}_3_1000K_li.txt",f"plots/{phase}_4_1200K_li.txt"]
'''
files_p = [f"plots/{phase}_0_400K_p.txt",f"plots/{phase}_1_500K_p.txt",f"plots/{phase}_2_600K_p.txt",f"plots/{phase}_3_700K_p.txt",f"plots/{phase}_4_800K_p.txt",f"plots/{phase}_5_900K_p.txt",f"plots/{phase}_6_1000K_p.txt"]
files_s = [f"plots/{phase}_0_400K_s.txt",f"plots/{phase}_1_500K_s.txt",f"plots/{phase}_2_600K_s.txt",f"plots/{phase}_3_700K_s.txt",f"plots/{phase}_4_800K_s.txt",f"plots/{phase}_5_900K_p.txt",f"plots/{phase}_6_1000K_p.txt"]
files_li = [f"plots/{phase}_0_400K_li.txt",f"plots/{phase}_1_500K_li.txt",f"plots/{phase}_2_600K_li.txt",f"plots/{phase}_3_700K_li.txt",f"plots/{phase}_4_800K_li.txt",f"plots/{phase}_5_900K_p.txt",f"plots/{phase}_6_1000K_p.txt"]
'''
files_li = files_li[0:2]+files_li[3:]
files_p = files_p[0:2]+files_p[3:]
files_s = files_s[0:2]+files_s[3:]
'''
T_colorbar = np.linspace(400,1000, 256)
Ts = list(range(400,1100,100))#[400,600,1000,1200]#[500,600,700, 800,900]
ratio = []
fig, ax = plt.subplots(figsize=(5,5))
for i,file in enumerate(files_p):
    data = np.loadtxt(files_p[i]).T
    data_s = np.loadtxt(files_s[i]).T
    data_li = np.loadtxt(files_li[i]).T
    ratio.append(data[-1,1]/data_li[-1,1])
    if i==0:
        ax.plot(data[::50, 0], data[::50, 1], label="P", color = cmap(i/4), linestyle="solid")
        ax.plot(data_s[::50, 0], data_s[::50, 1], label="S", color = cmap(i/4), linestyle="dotted")
    else:
        ax.plot(data[::50, 0], data[::50, 1], color=cmap(i/4), linestyle="solid")
        ax.plot(data_s[::50, 0], data_s[::50, 1], color=cmap(i/4), linestyle="dotted")
    plt.xlabel("t / ps")
    plt.ylabel("MSD / $\AA ^2$")
plt.legend()

norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')

plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap, ),
             ax=ax, orientation='horizontal',
             label="Temperature / K", ticks=[400,500,600,700,800, 900])
#plt.savefig("../../../Report/msd_ps.pdf")
plt.savefig(path+phase+"_new_msd_ps.png")
plt.figure()
for i,file in enumerate(files_p):
    data_li = np.loadtxt(files_li[i]).T
    plt.plot(data_li[:, 0], data_li[:, 1], label="T= "+str(Ts[i])+" K")
plt.legend()

print(ratio)
plt.savefig(path+phase+"_new_Na.png")
