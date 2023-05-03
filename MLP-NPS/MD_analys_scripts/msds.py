import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib import cm
timestep = 0.002
cmap = cm.get_cmap('coolwarm')
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.autolayout': True})

T_colorbar = np.linspace(400,800, 256)
folder = "/home/yli/database/triple_na/triple_na4p2s7/"
fig, ax = plt.subplots(figsize=(5, 5))
for i in range(0,6):
    data = np.loadtxt(folder+str(i)+"/msd.txt", ndmin=2)
    print(np.shape(data))
    plt.plot(data[:,0]*timestep, np.sum(data[:3000,1:], axis=1), color=cmap(i/4))
    plt.ylabel("Li msd / $\AA^2$")
    plt.xlabel("t / ps")
    ax.set_aspect(1. / ax.get_data_ratio())
norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')
'''
cmap = cm.get_cmap('coolwarm')
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.autolayout': True})

T_colorbar = np.linspace(400,800, 256)
folder = "li3ps4/plots/"
fig, ax = plt.subplots(figsize=(5, 5))
for i,T in enumerate([400, 500, 600]):
    data = np.loadtxt(folder+"beta_corr_lowT_"+str(T)+"_all_all.txt", ndmin=2)
    data = data.reshape((4,int(len(data)/4))).T
    print(np.shape(data))
    plt.plot(data[:3000,0], np.sum(data[:3000,1:], axis=1), color=cmap(i/4))
    plt.ylabel("Li msd / $\AA^2$")
    plt.xlabel("t / ps")
    ax.set_aspect(1. / ax.get_data_ratio())
norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')

plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap, ),
             ax=ax, orientation='horizontal',shrink=.8,pad = 0.2,
             label="Temperature / K", ticks=[400,500,600,700, 800])
plt.savefig("/home/huss/MA/Report/limsd_gamma.pdf")


fig, ax = plt.subplots(figsize=(5, 5))
for i,T in enumerate([400, 500, 600]):
    data = np.loadtxt(folder+"gamma_cr_"+str(T)+"_all_all.txt", ndmin=2)
    data = data.reshape((4,int(len(data)/4))).T
    print(np.shape(data))
    plt.plot(data[:3000,0], np.sum(data[:3000,1:], axis=1), color=cmap(i/4))
    plt.ylabel("Li msd / $\AA^2$")
    plt.xlabel("t / ps")
ax.set_aspect(1. / ax.get_data_ratio())
norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')

plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax, orientation='horizontal',shrink=.8,pad = 0.2,
             label="Temperature / K", ticks=[400,500,600,700, 800])
plt.savefig("/home/huss/MA/Report/limsd_beta.pdf")

fig, ax = plt.subplots(figsize=(5, 5))
for i,T in enumerate([500, 600, 700, 800]):
    data = np.loadtxt(folder+"alpha_true_"+str(T)+"_all_all.txt", ndmin=2)
    data = data.reshape((4, int(len(data)/4))).T
    print(np.shape(data))
    plt.plot(data[:3000,0], np.sum(data[:3000,1:], axis=1), color=cmap((i+1)/4))
    plt.ylabel("Li msd / $\AA^2$")
    plt.xlabel("t / ps")
#ax.set_aspect(1. / ax.get_data_ratio())
norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')

plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap, ),
             ax=ax, orientation='horizontal',shrink=.8, pad = 0.2,
             label="Temperature / K", ticks=[400,500,600,700, 800])
plt.savefig("/home/huss/MA/Report/limsd_alpha.pdf")


fig, ax = plt.subplots(figsize=(5, 5))
for i,T in enumerate([400, 500, 600, 700]):
    data = np.loadtxt(folder+"li7p3s11_"+str(T)+"_all_all.txt", ndmin=2)
    data = data.reshape((4, int(len(data)/4))).T
    print(np.shape(data))
    plt.plot(data[:3000,0], np.sum(data[:3000,1:], axis=1), color=cmap(i/4))
    plt.ylabel("Li msd / $\AA^2$")
    plt.xlabel("t / ps")
ax.set_aspect(1. / ax.get_data_ratio())
norm = mpl.colors.BoundaryNorm(T_colorbar, 257, extend='both')

plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap, ),
             ax=ax, orientation='horizontal',shrink=.8, pad = 0.2,
             label="Temperature / K", ticks=[400,500,600,700, 800])
plt.savefig("/home/huss/MA/Report/limsd_li7p3s11.pdf")

plt.show()
'''