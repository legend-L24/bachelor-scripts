#### script to get Li isodensity from trajectory
from ase import Atoms
from ase.io import read
from fourier import get_3d_density
import numpy as np
import logging


def remove_edgeatoms(atoms):
    """remove the atoms being very close to the boundary to avoid array mismatch problems"""
    atoms2 = Atoms(cell=atoms.cell, pbc=[True,True,True])
    for atom in atoms:
        if not any(np.abs(atom.position[i] - atoms.cell[i,i]) < 0.001 for i in range(3)):
            atoms2.append(atom)
    return atoms2


def map_traj_to_atoms(traj,element,freq,com=False,comrefelements=None):
    """generate atoms object containing all atoms of type el in traj"""
    ref = traj[0]
    ref.wrap()
    if com:
        refelements = comrefelements
        refatoms = [atom.index for atom in ref if atom.symbol in refelements]
    densatoms = [atom.index for atom in ref if atom.symbol==element]
    densatom = Atoms()

    for count,frame in enumerate(traj):
        if count%freq == 0:
            frame.wrap()
            if com:
                ref_com = ref[refatoms].get_center_of_mass()
                curr_com = frame[refatoms].get_center_of_mass()
                diff = ref_com - curr_com
                frame.translate(diff)
            for atom in frame[densatoms]:
                densatom.append(atom)
    densatom.wrap()
    densatom.set_cell(ref.get_cell())
    densatom.pbc=[True,True,True]
    densatom.wrap()
    densatom = remove_edgeatoms(densatom)
    return densatom


def sum_isovolume(dens_3d, th):
    """sum up the isovolume above a specific threshold"""
    shape=np.shape(dens_3d)
    return np.sum(dens_3d >= th)/(shape[0]*shape[1]*shape[2])


def get_vol_thresh_curve(full_traj, outfile, spacing=200, npoints=20):
    """default spacing: 200 fs * 200 = 40 ps"""
    traj = full_traj[::spacing]
    for atoms in traj:
        atoms.pbc = 1, 1, 1
        atoms.wrap()
    mapped_atoms = map_traj_to_atoms(traj, element='Li', freq=1, com=True, comrefelements=["P", "S"])
    dens = get_3d_density(mapped_atoms, elements="Li")
    dens = dens / np.linalg.norm(dens)
    threshs = np.linspace(min(dens.flatten()), max(dens.flatten()), npoints)
    vols = np.zeros(len(threshs))
    for j in range(len(threshs)):
        vols[j] = sum_isovolume(dens, threshs[j])
        with open(outfile, "a") as f:
            f.write("{thresh} {vol}\n".format(thresh=threshs[j], vol=vols[j]))


if __name__=="__main__":
    log = logging.getLogger("loop")
    for struct in ["Li3PS4", "Li7P3S11", "Li4P2S7"]:
        for T in [400, 500, 600, 700]:
            log.warning("T ={}".format(T))
            full_traj = read(r"E:\MA_backup\tabea_ma_data\md_runs\production_runs\statistics_MD\{}\{}\2\geom.xyz".format(struct, T), ":")
            get_vol_thresh_curve(full_traj, "vol_data/vol3_{}_{}.txt".format(struct, T))
