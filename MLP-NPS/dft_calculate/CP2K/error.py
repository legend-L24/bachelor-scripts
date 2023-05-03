#!/usr/bin/env python
from ase.io import read, write
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_squared_error


work_dir = "."
atoms = read(work_dir+ "/quip_train.xyz", "::")
atoms_= read(work_dir+ "/train.xyz", "::")
frc_dft = np.array([])
frc_gap = np.array([])
ener_dft = []
ener_gap = []
for idx, atom in enumerate(atoms[:-3:]):
    frc_dft = np.append(frc_dft, atom.get_forces())
    frc_gap = np.append(frc_gap, atom.get_array('force'))
    atom_len = len(atom.get_chemical_symbols())
    ener_dft.append(atoms_[idx].info['energy']/atom_len)
    ener_gap.append(atom.get_potential_energy()/atom_len)
ener_dft = np.array(ener_dft)
ener_gap = np.array(ener_gap)

print(np.min(ener_dft))
print(np.max(ener_dft))

lw2 =2
fs = 22
plt.figure(figsize=[18, 8], dpi=100)

ax = plt.subplot(121)
plt.title(r"(a) Potential Energy", fontsize = fs)
ax.spines['bottom'].set_linewidth(lw2)
ax.spines['left'].set_linewidth(lw2)
ax.spines['right'].set_linewidth(lw2)
ax.spines['top'].set_linewidth(lw2)
plt.plot([-20000,20000], [-20000,20000], "k--")
min_val = -7099.15
scale = 0.4
rmse = np.sqrt(mean_squared_error(ener_dft, ener_gap))
print(np.mean(ener_dft-ener_gap))
std = np.sqrt(np.var((ener_dft-ener_gap)**2))
plt.text(min_val +0.015/0.045*scale, min_val + 0.005/0.045*scale, "RMSE:\n" + r" %4.2e $\pm$ %4.2e" %(rmse, std), fontsize = fs)

plt.xticks(np.arange(min_val, min_val + scale+0.001, scale/2), fontsize = fs - 2)
plt.yticks(np.arange(min_val, min_val + scale+0.001, scale/2), fontsize = fs - 2)

plt.xlim(min_val, min_val+ scale)
plt.ylim(min_val, min_val+ scale)
plt.plot(ener_dft, ener_gap, "r.")
plt.ylabel(r"$\rm \Delta E_{DFT}$ (eV/atom)", fontsize = fs)
plt.xlabel(r"$\rm \Delta E_{GAP} $ (eV/atom)", fontsize = fs)

ax = plt.subplot(122)
plt.title(r"(b) Atomic Forces", fontsize = fs)
ax.spines['bottom'].set_linewidth(lw2)
ax.spines['left'].set_linewidth(lw2)
ax.spines['right'].set_linewidth(lw2)
ax.spines['top'].set_linewidth(lw2)

plt.plot(frc_gap, frc_dft, "b.")
rmse = np.sqrt(mean_squared_error(frc_dft, frc_gap))
scale = 10
print(np.mean(frc_dft-frc_gap))
std = np.sqrt(np.var((frc_dft-frc_gap)**2))
plt.text(-0.015/0.045*scale, -0.035/0.045*scale, "RMSE:\n" + r" %4.2e $\pm$ %4.2e" %(rmse, std), fontsize = fs)
plt.xticks(np.arange(scale*(-1), scale+0.001, scale), fontsize = fs - 2)
plt.yticks(np.arange(scale*(-1), scale+0.001, scale), fontsize = fs - 2)
plt.ylabel(r"$\rm\Delta F^i_{DFT}$ (eV/$\rm \AA$)", fontsize = fs)
plt.xlabel(r"$\rm\Delta F^i_{GAP}}$ (eV/$\rm \AA$)", fontsize = fs)

plt.plot([-20,20], [-20,20], "k--")
plt.xlim(-scale, scale)
plt.ylim(-scale, scale)

plt.savefig("energy.png",bbox_inches = 'tight')
