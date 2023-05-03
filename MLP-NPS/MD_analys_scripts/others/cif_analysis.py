from ase.io import read, write
import os
import matplotlib.pyplot as plt

vol_ls = []
vol_dic = dict({})
path = "/home/yli/test_small"
cif_ls = os.listdir(path)
for file in cif_ls:
    if file[-3:]!="cif": continue
    atom = read(os.path.join(path, file))
    vol_dic[file] = atom.get_volume()
    vol_ls.append(atom.get_volume())
vol_ls.sort()
name_ls = sorted(vol_dic)
picked_name = name_ls[50:100]
picked_name = [i for i in  picked_name if "2B" in i]
print(picked_name)
print(len(picked_name))

plt.plot(vol_ls)
plt.show()
