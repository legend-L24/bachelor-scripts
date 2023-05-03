# Poisson disc sampling in arbitrary dimensions via Bridson algorithm
# Implementation by Pavel Zun, pavel.zun@gmail.com
# BSD licence - https://github.com/diregoblin/poisson_disc_sampling

# -----------------------------------------------------------------------------
# Based on 2D sampling by Nicolas P. Rougier - https://github.com/rougier/numpy-book
# -----------------------------------------------------------------------------

import numpy as np
from scipy.special import gammainc
from copy import deepcopy
from ase import Atoms
from ase.visualize import view


# ToDo: do nice visualisation (use ase with x elements for visualization)
# ToDo: make periodic (wrap cells inside when close to boundary)

# Uniform sampling in a hyperspere
# Based on Matlab implementation by Roger Stafford
# Can be optimized for Bridson algorithm by excluding all points within the r/2 sphere


def hypersphere_volume_sample(center, radius, k=1):
    ndim = center.size
    x = np.random.normal(size=(k, ndim))
    ssq = np.sum(x ** 2, axis=1)
    fr = radius * gammainc(ndim / 2, ssq / 2) ** (1 / ndim) / np.sqrt(ssq)
    frtiled = np.tile(fr.reshape(k, 1), (1, ndim))
    p = center + np.multiply(x, frtiled)
    return p


# Uniform sampling on the sphere's surface
def hypersphere_surface_sample(center, radius, k=1):
    ndim = center.size
    vec = np.random.standard_normal(size=(k, ndim))
    vec /= np.linalg.norm(vec, axis=1)[:, None]
    p = center + np.multiply(vec, radius)
    return p


def squared_distance(p0, p1):
    return np.sum(np.square(p0 - p1))


def Bridson_sampling(dims=np.array([1.0, 1.0]), radius=0.05, k=30, initial_points=[],
                     hypersphere_sample=hypersphere_volume_sample):
    # References: Fast Poisson Disk Sampling in Arbitrary Dimensions

    #             Robert Bridson, SIGGRAPH, 2007
    ndim = dims.size
    # size of the sphere from which the samples are drawn relative to the size of a disc (radius)
    sample_factor = 2
    if hypersphere_sample == hypersphere_volume_sample:
        sample_factor = 2
    # for the surface sampler, all new points are almost exactly 1 radius away from at least one existing sample
    # eps to avoid rejection
    if hypersphere_sample == hypersphere_surface_sample:
        eps = 0.001
        sample_factor = 1 + eps

    def in_limits(p):
        return np.all(np.zeros(ndim) <= p) and np.all(p < dims)

    # Check if there are samples closer than "squared_radius" to the candidate "p"

    def in_neighborhood(p, n=2):
        indices = (p / cellsize).astype(int)
        indmin = np.maximum(indices - n, np.zeros(ndim, dtype=int))
        indmax = np.minimum(indices + n + 1, gridsize)
        # Check if the center cell is empty
        if not np.isnan(P[tuple(indices)][0]):
            return True
        a = []
        for i in range(ndim):
            a.append(slice(indmin[i], indmax[i]))
        if np.any(np.sum(np.square(p - P[tuple(a)]), axis=ndim) < squared_radius):
            return True

    def add_point(p):
        points.append(p)
        indices = (p / cellsize).astype(int)
        try:
            P[tuple(indices)] = p
        except:
            print("point not within bounds")

    def remove_initial(positions):
        new_pos = []
        for position in positions:
            append = True
            for old_position in initial_points:
                if np.linalg.norm(position - old_position) < 0.1: \
                        append = False
            if append:
                new_pos.append(position)
        return new_pos

    cellsize = radius / np.sqrt(ndim)
    gridsize = (np.ceil(dims / cellsize)).astype(int)
    # Squared radius because we'll compare squared distance
    squared_radius = radius * radius
    # Positions of cells
    P = np.empty(np.append(gridsize, ndim), dtype=np.float32)  # n-dim value for each grid cell
    # Initialise empty cells with NaNs
    P.fill(np.nan)
    points = []
    for initial_point in initial_points:
        add_point(initial_point)
    # add_point(np.random.uniform(np.zeros(ndim), dims))
    while len(points):
        i = np.random.randint(len(points))
        p = points[i]
        del points[i]
        Q = hypersphere_sample(np.array(p), radius * sample_factor, k)
        for q in Q:
            if in_limits(q) and not in_neighborhood(q):
                add_point(q)
    positions = P[~np.isnan(P).any(axis=ndim)]
    positions_nocr = remove_initial(positions)
    return positions_nocr


def Bridson_sampling_periodic(dims=np.array([1.0, 1.0]), radius=0.05, k=30, initial_points=[],
                              hypersphere_sample=hypersphere_volume_sample, axis=0):
    ndim = dims.size
    # size of the sphere from which the samples are drawn relative to the size of a disc (radius)
    sample_factor = 2
    if hypersphere_sample == hypersphere_volume_sample:
        sample_factor = 2
    # for the surface sampler, all new points are almost exactly 1 radius away from at least one existing sample
    # eps to avoid rejection
    if hypersphere_sample == hypersphere_surface_sample:
        eps = 0.001
        sample_factor = 1 + eps

    def get_periodic_dublicates(p):
        point_copies = [p]
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    if not (i==0 and j==0 and k==0):
                        p_ = deepcopy(p)
                        p_[0] += i*dims[0]
                        p_[1] += j*dims[1]
                        p_[2] += k*dims[2]
                        point_copies.append(p_)
        return point_copies

    def in_limits(p):
        return np.all(np.zeros(ndim) <= p) and np.all(p < dims)

    # Check if there are samples closer than "squared_radius" to the candidate "p"
    def in_neighborhood(p, n=2):
        plist = get_periodic_dublicates(p)
        for p in plist:
            indices = (p / cellsize).astype(int)
            indmin = np.maximum(indices - n, np.zeros(ndim, dtype=int))
            indmax = np.minimum(indices + n + 1, gridsize)
            # Check if the center cell is empty
            if not np.isnan(P[tuple(indices)][0]):
                return True
            a = []
            for i in range(ndim):
                a.append(slice(indmin[i], indmax[i]))
            if np.any(np.sum(np.square(p - P[tuple(a)]), axis=ndim) < squared_radius):
                return True

    def add_point(p):
        points.append(p)
        indices = (p / cellsize).astype(int)
        indices_copies = [indices]
        l = [0, 1, 2]
        l.pop(axis)
        point_copies = get_periodic_dublicates(p)
        for i in l:
            indices_ = deepcopy(indices)
            indices_[i] += initial_gridsize[i]
            indices_copies.append(indices_)
        indices_ = deepcopy(indices)
        indices_[l[0]] += initial_gridsize[l[0]]
        indices_[l[1]] += initial_gridsize[l[1]]
        indices_copies.append(indices_)
        for p_, indices_ in zip(point_copies, indices_copies):
            try:
                P[tuple(indices_)] = p_
            except:
                0

    def remove_initial(positions):
        new_pos = []
        for position in positions:
            append = True
            for old_position in initial_points:
                if np.linalg.norm(position - old_position) < 0.1: \
                        append = False
            if append:
                new_pos.append(position)
        return new_pos

    cellsize = radius / np.sqrt(ndim)
    gridsize = (np.ceil(dims / cellsize)).astype(int)
    initial_gridsize = deepcopy(gridsize)
    for i in range(3):
        gridsize[i] *= 2

    # Squared radius because we'll compare squared distance
    squared_radius = radius * radius
    # Positions of cells
    P = np.empty(np.append(gridsize, ndim), dtype=np.float32)  # n-dim value for each grid cell
    # Initialise empty cells with NaNs
    P.fill(np.nan)
    points = []
    for initial_point in initial_points:
        add_point(initial_point)
    if len(initial_points) == 0:
        add_point(np.array([3,3,3]))
    while len(points):
        i = np.random.randint(len(points))
        p = points[i]
        del points[i]
        Q = hypersphere_sample(np.array(p), radius * sample_factor, k)
        for q in Q:
            if in_limits(q) and not in_neighborhood(q):
                add_point(q)
    P = P[:initial_gridsize[0], :initial_gridsize[1], :initial_gridsize[2]]
    positions = P[~np.isnan(P).any(axis=ndim)]
    positions_nocr = remove_initial(positions)
    return positions_nocr
