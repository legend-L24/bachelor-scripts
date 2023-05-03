#!/home/ytli/anaconda3/envs/quip/bin/python
from ase.io import read, write
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.cluster import KMeans

ener_region = 2
N_cluster = 1



work_dir = "/home/yli"
atoms = read(work_dir+ "/quip_train.xyz", "::")
atoms_= read(work_dir+ "/test.xyz", "::")
frc_dft = np.array([])
frc_gap = np.array([])
ener_dft = []
ener_gap = []
'''
dist_ls = []

name_ls = ['P', 'P']

for idx, atom in enumerate(atoms[:]):
    print(atom[0].symbol,atom[1].symbol)
    if atom[0].symbol==name_ls[0] and atom[1].symbol==name_ls[1]:
        ener_dft.append(atoms_[idx].info['energy'])
        ener_gap.append(atom.get_potential_energy())
        dist_ls.append(atom.get_positions()[-1][0])
    elso: print("skip")

plt.title(name_ls[0]+'-'+name_ls[1]+" disorder curve")
plt.scatter(dist_ls, ener_gap, label='gap')
plt.scatter(dist_ls, ener_dft, label='dft')
plt.legend()
plt.savefig(name_ls[0]+'-'+name_ls[1]+" disorder curve"+".png")
plt.show()
'''
for idx, atom in enumerate(atoms[:]):
    frc_dft = np.append(frc_dft, atoms_[idx].get_array('forces'))
    frc_gap = np.append(frc_gap, atom.get_array('force'))
    atom_len = len(atom.get_chemical_symbols())
    ener_dft.append(atoms_[idx].info['energy']/atom_len)
    ener_gap.append(atom.get_potential_energy()/atom_len)

#ener_dft = [i%ener_region for i in ener_dft]
#ener_gap = [i%ener_region for i in ener_gap]

ener_dft = np.array(ener_dft)
ener_gap = np.array(ener_gap)
print("this is the whole number", ener_gap.shape)
print(np.min(ener_dft))
print(np.max(ener_dft))
print(np.min(ener_gap))
print(np.max(ener_gap))

n_cluster = N_cluster
kmeans_kwargs = {
    "init": "random",
    "n_init": 10,
    "max_iter": 300,
    "random_state": 42,
    "n_clusters": n_cluster
}

kmeans = KMeans(**kmeans_kwargs)
kmeans.fit(ener_dft.reshape(-1,1))
print(kmeans.labels_, kmeans.cluster_centers_)

ener_dft_new = []
ener_gap_new = []
for i, j in enumerate(kmeans.labels_):
    ener_dft_new.append(ener_dft[i]-kmeans.cluster_centers_[j])
    ener_gap_new.append(ener_gap[i]-kmeans.cluster_centers_[j])

ener_dft = np.array(ener_dft_new)
ener_gap = np.array(ener_gap_new)


lw2 = 2
fs = 32
plt.figure(figsize=[18, 8], dpi=100)

ax = plt.subplot(121)

#plt.title(r"Energy", fontsize = fs)
ax.spines['bottom'].set_linewidth(lw2)
ax.spines['left'].set_linewidth(lw2)
ax.spines['right'].set_linewidth(lw2)
ax.spines['top'].set_linewidth(lw2)
plt.plot([-20000,20000], [-20000,20000], "k--")
min_val = np.min(ener_dft)
min_val = np.min(ener_gap)
scale = (np.max(ener_dft)-np.min(ener_dft))*1.1
scale = (np.max(ener_gap)-np.min(ener_gap))*1.1
ener_dft = ener_dft - min_val
ener_gap = ener_gap - min_val
#scale = (np.max(ener_gap)-np.min(ener_gap))*2
#scale = 0.15
#scale = 0.4
scale*=1.2
rmse = np.sqrt(mean_squared_error(ener_dft, ener_gap))
print(np.mean(ener_dft-ener_gap))
std = np.sqrt(np.var((ener_dft-ener_gap)**2))*1.5
#plt.text(min_val +0.015/0.045*scale, min_val + 0.005/0.045*scale, "RMSE:\n" + r" %4.2e $\pm$ %4.2e" %(rmse, std), fontsize = fs-5)
#plt.text(min_val +0.015/0.045*scale, min_val + 0.005/0.045*scale, "RMSE:\n" + r" %4.2e (eV/atom)" %(rmse), fontsize = fs-5)
#plt.text(0.015/0.045*scale, 0.005/0.045*scale, "RMSE:\n" + r" %4.2e (meV/Atom)" %(rmse), fontsize = fs-5)
plt.text(0.015/0.045*scale, 0.005/0.045*scale, "RMSE:\n" + r" %.2f meV/Atom" %(rmse*1000), fontsize = fs-5)
#plt.xticks(np.arange(min_val, min_val + scale+0.001, scale/2), fontsize = fs - 2)
#plt.yticks(np.arange(min_val, min_val + scale+0.001, scale/2), fontsize = fs - 2)
print("this is scale:",scale)
#scale+=0.15
if 1:
    plt.xticks(np.arange(0, scale+0.001, scale/2), fontsize = fs - 2)
    plt.yticks(np.arange(0, scale+0.001, scale/2), fontsize = fs - 2)
if 0:
    plt.xticks(np.arange(0, scale+0.001, 0.01), fontsize = fs - 2)
    plt.yticks(np.arange(0, scale+0.001, 0.01), fontsize = fs - 2)
plt.xlim(-scale*0.05, scale)
plt.ylim(-scale*0.05, scale)
#plt.xlim(min_val-scale, min_val+ scale)
#plt.ylim(min_val-scale, min_val+ scale)
plt.plot(ener_dft, ener_gap, "r.",markersize=15)
plt.ylabel(r"$\rm \Delta E_{DFT}$ (meV/Atom)", fontsize = fs)
plt.xlabel(r"$\rm \Delta E_{GAP} $ (meV/Atom)", fontsize = fs)

ax = plt.subplot(122)
#plt.title(r"Forces", fontsize = fs)
ax.spines['bottom'].set_linewidth(lw2)
ax.spines['left'].set_linewidth(lw2)
ax.spines['right'].set_linewidth(lw2)
ax.spines['top'].set_linewidth(lw2)
plt.plot(frc_gap, frc_dft, "b.")
rmse = np.sqrt(mean_squared_error(frc_dft, frc_gap))
scale = 5
print(np.mean(frc_dft-frc_gap))
std = np.sqrt(np.var((frc_dft-frc_gap)**2))
#plt.text(-0.015/0.045*scale, -0.035/0.045*scale, "RMSE:\n" + r" %4.2e $\pm$ %4.2e" %(rmse, std), fontsize = fs-5)
plt.text(-0.015/0.045*scale, -0.035/0.045*scale, "RMSE:\n" + r" %4.2e eV/$\rm \AA$" %(rmse), fontsize = fs-5)
plt.xticks(np.arange(scale*(-1), scale+0.001, scale), fontsize = fs - 2)
plt.yticks(np.arange(scale*(-1), scale+0.001, scale), fontsize = fs - 2)
plt.ylabel(r"$\rm\Delta F^i_{DFT}$ (eV/$\rm \AA)$", fontsize = fs)
plt.xlabel(r"$\rm\Delta F^i_{GAP}}$ (eV/$\rm \AA)$", fontsize = fs)


plt.plot([-20,20], [-20,20], "k--")
plt.xlim(-scale, scale)
plt.ylim(-scale, scale)
#plt.subplots_adjust(wspace=0.5)
#plt.savefig("energy_lmp.png",bbox_inches = 'tight')
plt.tight_layout()
plt.show()

