#!/u/ytli/anaconda3/envs/quip/bin/python
import numpy as np
import matplotlib.pyplot as plt
path = "./rdf.txt"
#path = "./aver.rdf"
x_ls = []
ls_0 = []
ls_1 = []
ls_2 = []
for i in np.loadtxt(path, skiprows=1, comments="#", dtype=float):
    x_ls.append(i[1])
    ls_0.append(i[2])
    ls_1.append(i[6])
    ls_2.append(i[4])
plt.title("radius distribution function")
plt.xlim(min(x_ls), max(x_ls))
plt.plot(x_ls, ls_0, color="black", label="P-P")
plt.plot(x_ls, ls_1, color="red", label="S-S")
plt.plot(x_ls, ls_2, color="yellow", label="P-S")
plt.legend()
plt.xlabel("radius/A")
plt.ylabel("G(r)")
plt.savefig("rdf.png")
