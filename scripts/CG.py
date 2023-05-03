#!/usr/bin/python
import numpy as np
from scipy.special import erf
from warnings import warn
import itertools as itt
import pdb

class Grid():

    def __init__(self, start, stop, num=50, PBC = False):
        '''Create a grid with the dimensions specified in (start, stop) and with spacing num'''
        self.isCyclic = PBC
        self.buffer = None
        self.dim = len(start)
        if type(num) == type(0): num = (num,) * self.dim    #if num is integer, make it the same dim as start
        self.x = []
        self.dx = []
        for k in range(self.dim):
            x, dx = np.linspace(start[k], stop[k], num[k], retstep=True)
            self.x.append(x)
            self.dx.append(dx)
        self.dx = np.array(self.dx)
        self.shape = tuple([len(x) for x in self.x])
        self.box = np.array(start, dtype=np.float), np.array(stop, dtype=np.float)

    def cyclicShape(self, *stencils):
        n = np.zeros(self.dim)
        for s in stencils:
            n = np.fmax(n, s.n0)    # we assume left/right symmetric stencils
        n = n.astype(np.int) + 1    # "+1" to avoid indexing errors on PBC boundary
        start = np.array([x[0]  for x in self.x]) - n*self.dx
        stop  = np.array([x[-1] for x in self.x]) + n*self.dx
        newg = Grid(start, stop, num = np.array(self.shape) + 2*n)
        newg.PBCgrid = self
        self.cutoffDimension = n
        newg.cutoffDimension = n
        self.isCyclic = True
        for s in stencils: s.grid = newg
        return newg

    def foldBack(self, cube, pbcDims = None, inplace = False, dtype = None):
        if self.isCyclic: return
        if pbcDims is None: pbcDims = range(self.dim)
        if dtype is None: dtype = cube.dtype
        if dtype != cube.dtype: inplace = False
        pbc = self.PBCgrid
        n = self.cutoffDimension
        if max(abs(n - pbc.cutoffDimension)) > 0: warn('non-matching grids')
        slicemap = [ [ ( slice(None, n[k]),    slice(-n[k], None)  ),
                       ( slice(n[k], -n[k]),   slice(None)         ),
                       ( slice(-n[k], None),   slice(None, n[k])   ),
                     ] for k in pbcDims ]
        if inplace:
            cpbc = cube[[slice(k, -k) for k in n]]
            center = (slice(None),) * self.dim
        else:
            if self.buffer is None:
                self.buffer = np.zeros(pbc.shape, dtype = dtype)
            else:
                self.buffer.fill(0)
            cpbc = self.buffer
        for m in itt.product(*slicemap):
            src, dest = zip(*m)
            if inplace and dest == center: continue # exclude non-wrapped part
#           print src, ' -> ', dest
            cpbc[dest] += cube[src]
        return pbc, cpbc
                  
class Stencil():

    def __init__(self, grid, norm = 1):
        self.grid = grid
        self.norm = float(norm)

    def Gaussian(self, alpha, epsilon = 1.0e-5):
        grid = self.grid
        a = np.sqrt(abs(alpha))
        width = np.sqrt(-np.log(epsilon**(1./grid.dim)/grid.dx)/alpha)
        nd = np.ceil(width/grid.dx).astype(np.int)
        stencil = np.ones(2*nd+1)
        for k in range(grid.dim):
            n = nd[k]
            dx = grid.dx[k]
            x = np.linspace(-n*dx, n*dx, num=2*n+1)
            f = erf(a * (x + dx/2.)) - erf(a * (x - dx/2.))
            idx = [None,]*grid.dim
            idx[k] = slice(None)
            stencil = np.multiply(stencil, f[tuple(idx)], stencil) #Changed due to FutureWarning of numpy
            #stencil = np.multiply(stencil, f[idx], stencil)
        self.n0 = nd    # distance of center w.r.t. "lower-left" corner
        self.shape = stencil.shape
        self.stencil = stencil
        self.normalize()
        return a, width, nd, stencil

    def normalize(self, norm = None):
        if norm is not None: self.norm = float(norm)
        stencil = self.stencil
        norm = self.norm/np.sum(stencil)
        self.stencil = np.multiply(stencil, norm, stencil)
        if max(self.shape) > max(self.grid.shape):
            warn('stencil is too large for grid')

    def move(self, r0):
        grid = self.grid
        idx = [np.searchsorted(grid.x[k], r0[k]) for k in range(grid.dim)]
        weight = np.array([grid.x[k][idx[k]] - r0[k] for k in range(grid.dim)])/grid.dx
        return np.array(idx), weight

    def corners(self, r0):
        idx, weight = self.move(r0)
        result = []
        for corner in itt.product(*zip(zip(idx-1, 1-weight), zip(idx, weight))):
            result.append( (np.array  ([p[0] for p in corner]),
                            np.product([p[1] for p in corner])) )
        return result

    def place(self, cube, centers, scale = 1.):
        scaled = np.zeros(self.stencil.shape)
        count = 0
        for r0 in centers:
            count += 1
            for corner, weight in self.corners(r0):
                start = corner - self.n0 + 1
                domain = [slice(*r) for r in zip(start, start + self.shape)]
                #print(count,cube[tuple(domain)].shape, tuple(domain))
                cube[tuple(domain)] += np.multiply(self.stencil, weight*scale, scaled) #Changed due to np FutureWarning
                #cube[domain] += np.multiply(self.stencil, weight*scale, scaled)

    def place_value(self, cube, centers, values = None):
        if values is None: values = np.ones(len(centers))
        scaled = np.zeros(self.stencil.shape)
        for r0,val in zip(centers, values):
            for corner, weight in self.corners(r0):
                start = corner - self.n0 + 1
                domain = [slice(*r) for r in zip(start, start + self.shape)]
                cube[tuple(domain)] += np.multiply(self.stencil, weight*val, scaled) #Changed due to np FutureWarning

                #cube[domain] += np.multiply(self.stencil, weight*val, scaled)



def test3D(num=50, pbcDims = None):
    a = 5
    dim = 3
    g = Grid((-a,)*dim, (a,)*dim, num=num)
    s1 = Stencil(g)
    s1.Gaussian(0.5)
    s2 = Stencil(g)
    s2.Gaussian(1.5)
    g = g.cyclicShape(s1, s2)
    c = np.zeros(g.shape)   # use c.fill(0) to clear the cube!
    s1.place(c,((1.2,1.0,-0.1), (5,5,5), (-5,4.9,0),))
    s2.place(c,((0,0,0),))
    g, c = g.foldBack(c, pbcDims = pbcDims, dtype=np.float32)
#   g, c = g.foldBack(c, pbcDims = pbcDims, inplace=True)
#   g, c = g.foldBack(c, pbcDims = pbcDims)
#   surf = mayaviContourSurface(g, c)
#   for i in range(nit):
#       pts, cells, norms = surf(c)
    return g, c

def testcorners(a = 5, dim = 3, num = 50, single = True, pbcDims = None):
    g = Grid((-a,)*dim, (a,)*dim, num=50)
    s = Stencil(g)
    s.Gaussian(0.5)
    g = g.cyclicShape(s)
    c = np.zeros(g.shape)   # use c.fill(0) to clear the cube!
    corners = itt.product(*[(-a,a),]*g.dim)
    if single: corners = itt.islice(corners, 1)  # only one corner for PBC check
    s.place(c, corners)
    return g.foldBack(c, pbcDims = pbcDims)

if False and __name__ == '__main__':
    import multiprocessing
    import subprocess

    pool = multiprocessing.Pool(None)
    tasks = range(96)
    results = []
    r = pool.map_async(performanceTest, tasks, callback=results.append)
    r.wait() # Wait on the results

    from MarchingCubes.CG import testcorners
    from MarchingCubes.graph import contourPlot
    g,c=testcorners(pbcDims=[1,])
    contourPlot(g,c)
    plt.show()
