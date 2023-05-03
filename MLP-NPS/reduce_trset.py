"""
This code aims to delte similar structures from training dataset based on SOAP descriptor
"""


from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
from dscribe.descriptors.soap import SOAP
from sklearn.manifold import TSNE
from sklearn.decomposition import KernelPCA
from ase.visualize import view

def get_soap_descriptors(atoms, global_soap=True):
    soap = SOAP(rcut=6, n_max=8, l_max=3, sigma=.5, species=[11, 15, 16], periodic=True)
    local_soap = soap.create(atoms)
    if global_soap:
        return np.average(local_soap, axis=0)
    else:
        return local_soap


def get_tsne(descriptorarray):
    return TSNE(learning_rate="auto").fit_transform(descriptorarray)


def get_kpca(descriptorarray):
    return KernelPCA(n_components=2, kernel="rbf").fit_transform(descriptorarray)


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
    return sample_inds


def get_new_pcas(pca, indexes):
    """only visualization"""
    new_pca = []
    for index in indexes:
        new_pca.append(pca[index])
    return np.array(new_pca)


def get_new_dataset(dataset, indexes):
    new_dataset = []
    for index in indexes:
        new_dataset.append(dataset[index])
    return(new_dataset)

def get_rest_dataset(dataset, indexes):
    rest_dataset = []
    for index in range(0, len(dataset)):
        if index in indexes: continue
        rest_dataset.append(dataset[index])
    return(rest_dataset)

# read dataset
trset_0 = read("/home/yli/train.xyz", ":")


len_ls = []
for atom in trset_0:
    len_ls.append(len(atom))
plt.scatter(range(0, len(len_ls)),len_ls)
plt.show()

#write("tabea_data.xyz", trset_0[:150])
#trset_0 = trset_0[150:]
#print("this is other: ",len(trset_0))
picked_ls = [i for i,j in enumerate(trset_0) if len(j)>3 and i>153]
other_ls = [i for i,j in enumerate(trset_0) if len(j)<3 or i<153]
print("this is dimers legenth: ", len(picked_ls),len(other_ls))
trset = []
for idx in picked_ls:
    trset.append(trset_0[idx])
extra_trset = []
for idx in other_ls:
    extra_trset.append(trset_0[idx])
write("extra_structure.xyz", extra_trset)
print("this is dimers legenth: ",len([i for i,j in enumerate(trset) if len(j)>3]))
print("this is real length", len(trset))
#trset = [i for i,j in enumerate(trset) if len(j)>3]
#trset = trset[[i for i,j in enumerate(trset) if len(j)>3]]
# compute global soap
all_soaps = []
for i in range(len(trset)):
    all_soaps.append(get_soap_descriptors(trset[i]))
all_soaps = np.array(all_soaps)
# reduce to two dimensions
pca = get_kpca(all_soaps)
plt.scatter(pca[:, 0], pca[:,1])
indexes = fps(pca, 410)
new_pca = get_new_pcas(pca, indexes)
plt.scatter(new_pca[:,0], new_pca[:,1], marker="x")
plt.show()
new_dataset = get_new_dataset(trset, indexes)
rest_dataset = get_rest_dataset(trset, indexes)
#view(new_dataset)
write("train_p2s6_again.xyz", new_dataset)
write("train_p2s6_again_delete.xyz", rest_dataset)
# do farthest point sampling and get indexes back

