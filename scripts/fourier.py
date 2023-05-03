import numpy as np
from ase.io import read
import CG_christoph.CG as cg 
from scipy import fftpack, ndimage


def get_3d_density(atoms,elements,res=None,alpha=None,numfact=None,num=None):
    cell  = atoms.cell
    start = [0.,0.,0.]
    stop  = [cell[i,i] for i in range(len(cell))]
    if res == None:
        res = 0.1
    else:
        res = res

    if alpha == None:
        alpha = 0.25
    else:
        alpha = alpha
    if num is None:
        if numfact == None:
            num   = ((np.floor(stop)-np.round(start))/res+1).tolist()
            num   = [int(n) for n in num]
        else:
            num   = ((np.floor(stop)-np.round(start))/res+1).tolist()
            num   = [int(n)*numfact for n in num]
    
    centers = [atom.position for atom in atoms if atom.symbol in elements]

    g = cg.Grid(start=start,stop=stop,num=num)
    s = cg.Stencil(g)
    s.Gaussian(alpha)
    g = g.cyclicShape(s)
    c = np.zeros(g.shape, dtype=np.float32)
    s.place(c, list(centers))
    g,c_el = g.foldBack(c, pbcDims=None, dtype=np.float32)

    return c_el

