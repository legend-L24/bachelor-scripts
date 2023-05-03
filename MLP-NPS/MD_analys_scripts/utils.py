from ase.io import read, write
import numpy as np
from ase import Atoms


def get_pbc_cube(geom):
    geom_all = geom.copy()
    factors = get_mult_factors()
    for factor in factors:
        if not factor == [0,0,0]:
            add_geom = geom.copy()
            add_geom.positions+=np.dot(geom.cell, np.array(factor))
            geom_all.extend(add_geom)
    return geom_all


def get_composition(geom):
    p2s6 = 0
    p2s7 = 0
    geom_p = select_element_atoms(geom, "P")
    geom_p_all = get_pbc_cube(geom_p)
    for i1 in range(len(geom_p)):
        for i2 in range(len(geom_p_all)):
            if i1 != i2:
                dist = np.linalg.norm(geom_p[i1].position - geom_p_all[i2].position)
                if dist < 2.5:
                    p2s6+=1
                    continue
                elif 2.5 < dist < 4.2:
                    p2s7+=1
                    continue
    return len(geom_p)-p2s6-p2s7, p2s7, p2s6


def select_element_atoms(atoms, element):
    """select all atoms of one type and returns a atoms object containing only those atoms"""
    atoms2 = Atoms()
    for i, symb in enumerate(atoms.symbols):
        if symb == element:
            atoms2.append(atoms[i])
    atoms2.cell = atoms.cell
    return atoms2


def get_mult_factors():
    """static function, mult factors for pbc"""
    factors = []
    for x in [-1,0,1]:
        for y in [-1,0,1]:
            for z in [-1,0,1]:
                if not [x,y,z] == [0,0,0]:
                    factors.append([x,y,z])
    return factors