import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean


def dist_from_point(ponto,lista):
    return [ euclidean(ponto,lista[j]) for j in range(len(lista)) ]


def get_max_dist_point(lista_ds):
    ds_max = max(lista_ds)
    idx = lista_ds.index(ds_max)
    return pts[idx]


def fps(pts, N, K):
    farthest_pts = [0] * K

    P0 = pts[np.random.randint(0, N)]
    farthest_pts[0] = P0
    ds0 = dist_from_point(P0, pts)

    ds_tmp = ds0
    for i in range(1, K):
        farthest_pts[i] = get_max_dist_point(ds_tmp)
        ds_tmp2 = dist_from_point(farthest_pts[i], pts)
        ds_tmp = [min(ds_tmp[j], ds_tmp2[j]) for j in range(len(ds_tmp))]
        # print ('P[%d]: %s' % (i, farthest_pts[i]))
    return farthest_pts


if __name__ == "__main__":
    N, K = 80, 40
    pts = np.random.random_sample((N, 2))
    plt.scatter(pts[:,0], pts[:,1])
    fp = np.array(op(pts, N, K))
    plt.scatter(fp[:,0], fp[:,1])
    plt.show()
