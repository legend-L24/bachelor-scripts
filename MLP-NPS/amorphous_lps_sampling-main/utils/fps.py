import numpy as np
from copy import deepcopy
from ase import Atoms
from ase.visualize import view


def fps(points, n_samples):
    """
    points: [N, 3] array containing the whole point cloud
    n_samples: samples you want in the sampled point cloud typically << N
    """
    points = np.array(points)
    points_left = np.arange(len(points))  # [P]
    sample_inds = np.zeros(n_samples, dtype='int')  # [S]
    dists = np.ones_like(points_left) * float('inf')  # [P]
    selected = 0
    sample_inds[0] = points_left[selected]
    points_left = np.delete(points_left, selected)  # [P - 1]
    for i in range(1, n_samples):
        last_added = sample_inds[i - 1]
        dist_to_last_added_point = (
                (points[last_added] - points[points_left]) ** 2).sum(-1)  # [P - i]
        dists[points_left] = np.minimum(dist_to_last_added_point,
                                        dists[points_left])  # [P - i]
        selected = np.argmax(dists[points_left])
        print(selected)
        sample_inds[i] = points_left[selected]
        points_left = np.delete(points_left, selected)
    return points[sample_inds], points[points_left]


def conditional_fps(points, n_samples, initial_inds, axis, dims):

    def get_periodic_dublicates(points):
        points_periodic = []
        for i, ifact in zip([0]*3, [-1, 0, 1]):
            for j, jfact in zip([1]*3, [-1, 0, 1]):
                for k, kfact in zip([2]*3, [-1, 0, 1]):
                    p_ = deepcopy(points)
                    p_[:, i] += ifact*dims[i]
                    p_[:, j] += jfact*dims[j]
                    p_[:, k] += kfact*dims[k]
                    points_periodic.append(p_)
        return points_periodic

    """
    points: [N, 3] array containing the whole point cloud
    n_samples: samples you want in the sampled point cloud typically << N
    initial_inds: initial indexes (in this case the crystal, to avoid that points at the edge of the crystal are attributed as P)
    """
    points = np.array(points)
    points_left = np.arange(len(points))  # [P]
    sample_inds = np.zeros(n_samples+len(initial_inds), dtype='int')  # [S]
    sample_inds[:len(initial_inds)] = initial_inds
    dists = np.ones_like(points_left) * float('inf')  # [P]
    # reverse indexes to be able to remove them from points_left with their index
    initial_inds.reverse()
    for initial_ind in initial_inds:
        points_left = np.delete(points_left, initial_ind)  # [P - 1]
    # add periodic images
    points_periodic = get_periodic_dublicates(points)
    for i in range(len(initial_inds), n_samples+len(initial_inds)):
        last_added = sample_inds[i - 1]
        dist_to_last_added_point_list = []
        for j in range(len(points_periodic)):
            dist_to_last_added_point_list.append((
                    (points[last_added] - points_periodic[j][points_left]) ** 2).sum(-1))  # [P - i]
        dist_to_last_added_point_list = np.array(dist_to_last_added_point_list)
        dist_to_last_added_point = np.amin(dist_to_last_added_point_list, axis=0)
        dists[points_left] = np.minimum(dist_to_last_added_point,
                                        dists[points_left])  # [P - i]
        selected = np.argmax(dists[points_left])
        sample_inds[i] = points_left[selected]
        points_left = np.delete(points_left, selected)
    return points[sample_inds[-n_samples:]], points[points_left]


if __name__=="__main__":
    print("test run")