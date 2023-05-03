# Script to get partial RDFs given an ASE atoms object

import numpy as np
import math
from ase.io import read, write
from ase.visualize import view


def get_partial_rdf(atoms, rmax, nbins, el1="P", el2="Li"):
    '''
    Parameters:
    -----------
    atoms:      ASE atoms object
    rmax:       float, maximum distance in \AA to perform rdf
    nbins:      int, number of bins
    center_idx: list, list of indices which to reference over for rdf
                (e.g. mask of all Al atoms)
    rest_idx:   list, list of indices which to account for in rdf
                (e.g. mask of all Li atoms)

    Returns:
    --------
    rdf:        np.array, with rdf values
    dists:      np.array, with distance values
    '''

    def get_el_atoms(atoms, el):
        idxlist = []
        for i, atom in enumerate(atoms):
            if atom.symbol==el:
                idxlist.append(i)
        return idxlist

    import numpy as np
    import math

    dm = atoms.get_all_distances(mic=True)
    rdf = np.zeros(nbins + 1)
    dr = float(rmax / nbins)

    center_idx = np.array(get_el_atoms(atoms, el1))
    rest_idx = np.array(get_el_atoms(atoms, el2))
    for i in center_idx:
        for j in np.array(rest_idx):
            rij = dm[i][j]
            if rij != 0:
                index = int(math.ceil(rij / dr))
                if index <= nbins:
                    rdf[index] += 1

    dists = []
    for i in range(1, nbins + 1):
        rrr = (i - 0.5) * dr
        dists.append(rrr)
        # Normalize with spherical approximation
        outer = 4 / 3 * np.pi * (rrr + dr) ** 3
        inner = 4 / 3 * np.pi * (rrr) ** 3
        vol = outer - inner
        rdf[i] /= vol

    return np.array(rdf[1:]), np.array(dists)